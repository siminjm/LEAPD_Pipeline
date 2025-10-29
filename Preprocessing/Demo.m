% ==========================================================
% EEG Preprocessing Demo
% Author: Simin Jamshidi (University of Iowa)
% ==========================================================

clear; close all; clc;

fprintf('\n=========================\n');
fprintf('   EEG Preprocessing Demo\n');
fprintf('=========================\n\n');

%% ----------------------------------------------------------
% Step 1 — Generate synthetic multi-component EEG
% ----------------------------------------------------------
fprintf('Step 1: Generating synthetic EEG signal...\n');
Fs = 500;                % Sampling rate (Hz)
T = 5;                   % Duration (s)
t = (0:1/Fs:T-1/Fs)';    % Time vector
nCh = 8;                 % Number of channels

% --- Simulated rhythms ---
alpha = sin(2*pi*10*t);                    % 10 Hz
theta = 0.5*sin(2*pi*6*t);                 % 6 Hz
beta  = 0.2*sin(2*pi*20*t);                % 20 Hz
baseEEG = alpha + theta + beta;

% --- Channel-specific variation ---
X = zeros(length(t), nCh);
for ch = 1:nCh
    phaseShift = rand*2*pi;
    X(:,ch) = baseEEG .* (0.8+0.4*rand) + 0.05*sin(2*pi*60*t+phaseShift); % +60 Hz noise
end

% --- Add blink artifact on channels 1–2 ---
blink = exp(-((t-2.5)/0.12).^2) * 5;  % strong transient
X(:,1:2) = X(:,1:2) + blink;

% --- Add small Gaussian noise ---
X = X + 0.1*randn(size(X));

labels = arrayfun(@(i)sprintf('Ch%d', i), 1:nCh, 'UniformOutput', false);
fprintf('Synthetic EEG generated: %d samples × %d channels\n\n', size(X,1), size(X,2));


%% ----------------------------------------------------------
% Step 2 — Run preprocessing pipeline
% ----------------------------------------------------------
fprintf('Step 2: Running preprocessing pipeline...\n');
[X_clean, labels_clean, report, savePath] = pipeline_preprocessing(X, Fs, labels);

%% ----------------------------------------------------------
% Step 3 — Visualization
% ----------------------------------------------------------
fprintf('Step 3: Plotting before/after comparison...\n');

chToPlot = 1;   % visualize channel 1
t_clean = (0:length(X_clean)-1)/Fs;
t_raw   = (0:length(X)-1)/Fs;

figure('Color','w','Position',[200 150 1200 700]);
sgtitle('EEG Preprocessing — Before vs After','FontSize',14,'FontWeight','bold');

% --- Subplot 1: Raw EEG ---
subplot(3,1,1);
plot(t_raw, X(:,chToPlot),'k'); hold on;
xline(2.5,'r--','Blink artifact','LabelOrientation','horizontal','LabelVerticalAlignment','bottom');
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
title(sprintf('Raw EEG (%s)', labels{chToPlot}));
ylim([-2 5]); xlim([0 T]);
grid on;

% --- Subplot 2: Cleaned EEG ---
subplot(3,1,2);
plot(t_clean, X_clean(:,chToPlot),'b');
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
title(sprintf('Cleaned EEG (%s)', labels_clean{min(chToPlot,numel(labels_clean))}));
ylim([-2 5]); xlim([0 T]);
grid on;

% --- Subplot 3: Power Spectrum (before vs after) ---
subplot(3,1,3);
NFFT = 2^nextpow2(size(X,1));
f = Fs/2*linspace(0,1,NFFT/2+1);

rawFFT = abs(fft(X(:,chToPlot),NFFT)/size(X,1));
cleanFFT = abs(fft(X_clean(:,chToPlot),NFFT)/size(X_clean,1));

plot(f, 2*rawFFT(1:NFFT/2+1), 'k', 'LineWidth',1); hold on;
plot(f, 2*cleanFFT(1:NFFT/2+1), 'b', 'LineWidth',1.2);
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title(sprintf('Power Spectrum (%s)', labels{chToPlot}));
xlim([0 100]);
legend({'Raw','Cleaned'});
grid on;

%% ----------------------------------------------------------
% Step 4 — Summary
% ----------------------------------------------------------
fprintf('\n--- Preprocessing Summary ---\n');
fprintf('Line noise removed: %d\n', report.line_noise_removed);

if isfield(report, 'artifact_ics')
    fprintf('Artifact ICs removed: %d\n', numel(report.artifact_ics));
elseif isfield(report, 'removed_ICs')
    fprintf('Artifact ICs removed: %d\n', numel(report.removed_ICs));
else
    fprintf('Artifact ICs removed: (field missing)\n');
end

if isempty(report.bad_channels)
    fprintf('No noisy channels detected.\n');
else
    fprintf('Noisy channels removed: %s\n', strjoin(report.bad_channels, ', '));
end

if exist('savePath','var') && isfile(savePath)
    fprintf('Cleaned EEG saved at:\n  %s\n', savePath);
end

fprintf('-----------------------------\n');
fprintf('Demo complete.\n\n');
