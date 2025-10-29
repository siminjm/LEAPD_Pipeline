function [X_clean, labels_clean, report, savePath] = pipeline_preprocessing(X, Fs, labels)
% ===============================================================
% EEG Preprocessing Pipeline (Final Stable ICLabel Version)
% Author: Simin Jamshidi (University of Iowa)
% ===============================================================

fprintf('\n==========================================\n');
fprintf('        EEG Preprocessing Pipeline\n');
fprintf('==========================================\n');

% ---------- Configuration ----------
cfg.remove_line_noise = false;     % true if raw EEG, false if LEAPD already notched 60 Hz
cfg.detect_noisy_channels = true;
cfg.use_ica = true;
cfg.save_results = true;
cfg.save_dir = fullfile(pwd, 'cleaned_data');
if ~exist(cfg.save_dir, 'dir'), mkdir(cfg.save_dir); end

report = struct();

% --------------------------------------------------------------
% Step 1 — Line-noise removal (optional)
% --------------------------------------------------------------
if cfg.remove_line_noise
    fprintf('--- Step 1: Removing line noise ---\n');
    X_denoised = remove_line_noise(X, Fs);
    report.line_noise_removed = 1;
else
    fprintf('--- Step 1: Skipping line-noise removal (already handled in LEAPD) ---\n');
    X_denoised = X;
    report.line_noise_removed = 0;
end

% --------------------------------------------------------------
% Step 2 — Noisy-channel detection
% --------------------------------------------------------------
fprintf('--- Step 2: Detecting noisy channels ---\n');
if cfg.detect_noisy_channels
    [badChans, chStd] = detect_noisy_channels(X_denoised, labels);
    if ~isempty(badChans)
        fprintf('Detected noisy channels: %s\n', strjoin(labels(badChans), ', '));
        X_noisyRemoved = X_denoised;
        X_noisyRemoved(:, badChans) = [];
        labels_clean = labels;
        labels_clean(badChans) = [];
    else
        fprintf('No noisy channels detected.\n');
        X_noisyRemoved = X_denoised;
        labels_clean = labels;
    end
    report.bad_channels = labels(badChans);
else
    X_noisyRemoved = X_denoised;
    labels_clean = labels;
    report.bad_channels = {};
end

% --------------------------------------------------------------
% Step 3 — ICA + ICLabel artifact removal
% --------------------------------------------------------------
if cfg.use_ica
    fprintf('--- Step 3: Running ICA + ICLabel artifact removal ---\n');

    % Add EEGLAB path if needed
    eeglab_path = fullfile(pwd, 'eeglab2025.1.0');
    if exist(eeglab_path, 'dir'), addpath(eeglab_path); end
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui'); %#ok<ASGLU>

    % ----------------------------------------------------------
    % Import data (EEGLAB expects samples x channels)
    % ----------------------------------------------------------
    EEG = pop_importdata('dataformat','array', ...
                         'nbchan', numel(labels_clean), ...
                         'data', X_noisyRemoved', ... % ensure proper orientation
                         'srate', Fs);

    EEG.nbchan = numel(labels_clean);
    EEG.data = double(EEG.data);
    EEG.pnts = size(EEG.data, 2);
    EEG.trials = 1;
    EEG.xmin = 0;
    EEG.xmax = (EEG.pnts-1)/EEG.srate;

    % ----------------------------------------------------------
    % Assign valid synthetic channel locations
    % ----------------------------------------------------------
    nch = numel(labels_clean);
    tmpChans = repmat(struct( ...
        'labels','', 'X',[], 'Y',[], 'Z',[], ...
        'theta',[], 'radius',[], 'sph_theta',[], 'sph_phi',[], ...
        'sph_radius',[], 'type','EEG', 'ref','common', 'urchan',[]), 1, nch);

    for i = 1:nch
        tmpChans(i).labels = labels_clean{i};
        tmpChans(i).theta  = (i-1)*(360/nch);
        tmpChans(i).radius = 0.5;
        tmpChans(i).X = cosd(tmpChans(i).theta);
        tmpChans(i).Y = sind(tmpChans(i).theta);
        tmpChans(i).Z = 0;
        tmpChans(i).sph_theta  = tmpChans(i).theta;
        tmpChans(i).sph_phi    = 0;
        tmpChans(i).sph_radius = 1;
        tmpChans(i).urchan     = i;
    end
    EEG.chanlocs = tmpChans;

    % ----------------------------------------------------------
    % Minimal but valid chaninfo (prevents insertchans crash)
    % ----------------------------------------------------------
    EEG.chaninfo = struct();
    EEG.chaninfo.filename   = '';
    EEG.chaninfo.plotrad    = 0.5;
    EEG.chaninfo.shrink     = 1;
    EEG.chaninfo.nosedir    = '+X';
    EEG.chaninfo.ref        = 'common';
    EEG.chaninfo.chantype   = {{'EEG'}};
    EEG.chaninfo.nodatchans = [];
    EEG.chaninfo.nodatchanscom = '';

    % ----------------------------------------------------------
    % Finalize EEG structure and register geometry
    % ----------------------------------------------------------
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = pop_chanedit(EEG, 'eval','chans = EEG.chanlocs;');
    EEG = eeg_checkset(EEG);

    % ----------------------------------------------------------
    % Run ICA and ICLabel
    % ----------------------------------------------------------
    EEG = pop_runica(EEG, 'extended',1,'interrupt','off');
    EEG = eeg_checkset(EEG);

    EEG = iclabel(EEG);

    % ----------------------------------------------------------
    % Remove artifact ICs automatically
    % ----------------------------------------------------------
    badICs = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) < 0.5 & ...
        (EEG.etc.ic_classification.ICLabel.classifications(:,2) > 0.3 | ...
         EEG.etc.ic_classification.ICLabel.classifications(:,3) > 0.3 | ...
         EEG.etc.ic_classification.ICLabel.classifications(:,4) > 0.3));

    fprintf('Removing %d artifact ICs (blink/muscle/heart/etc.)...\n', numel(badICs));
    if ~isempty(badICs)
        EEG = pop_subcomp(EEG, badICs, 0);
    end
    EEG = eeg_checkset(EEG);

    X_clean = EEG.data';
    report.artifact_ics = badICs;
else
    fprintf('Skipping ICA + ICLabel step.\n');
    X_clean = X_noisyRemoved;
    report.artifact_ics = [];
end

% --------------------------------------------------------------
% Step 4 — Save cleaned data
% --------------------------------------------------------------
if cfg.save_results
    tstamp = datestr(now, 'yyyymmdd_HHMMSS');
    savePath = fullfile(cfg.save_dir, ['EEG_cleaned_' tstamp '.mat']);
    save(savePath, 'X_clean', 'labels_clean', 'report', 'Fs');
    fprintf('\n Cleaned EEG saved: %s\n', savePath);
else
    savePath = '';
end

% --------------------------------------------------------------
% Summary
% --------------------------------------------------------------
fprintf('\n--- Preprocessing Summary ---\n');
fprintf('Line noise removed: %d\n', report.line_noise_removed);
fprintf('Artifact ICs removed: %d\n', numel(report.artifact_ics));
if isempty(report.bad_channels)
    fprintf('No noisy channels detected.\n');
else
    fprintf('Noisy channels removed: %s\n', strjoin(report.bad_channels, ', '));
end
fprintf('Cleaned EEG saved at:\n  %s\n', savePath);
fprintf('-----------------------------\n');

end
