function results = main_train(cfg)
% MAIN_TRAIN - Unified training for classification and correlation.

% cfg fields:
%   .mode            : "classification" | "correlation"
%   .data_train      : path to training EEG .mat (struct or nested cells)
%   .labels_file     : table (xlsx/csv) with target variable for correlation
%   .channel_list    : cellstr of channel names to consider (optional)
%   .exclude_chans   : cellstr of noisy channels to skip (e.g., {'FT9','FT10','TP9','TP10'})
%   .f1_grid         : vector (e.g., 0.1:0.1:0.5)
%   .f2_grid         : vector (e.g., 1:1:100)
%   .orders          : vector (e.g., 2:10)
%   .save_dir        : 'results/train_results'
%   .Fs              : sampling rate (default 500)
%   .target_group    : 'group1'|'group2' (for correlation; target to correlate)
%   .kfold           : 'loocv' or integer K (for classification)
%   .is_norm_proj    : 0 or 1 (projection mode)
%   .is_preprocessed : true/false - indicate if data is already preprocessed
%   .detect_additional_noisy : true/false - auto-detect noisy channels

arguments
    cfg.mode (1,1) string {mustBeMember(cfg.mode,["classification","correlation"])}
    cfg.data_train (1,1) string
    cfg.labels_file (1,1) string = ""
    cfg.channel_list = []
    cfg.exclude_chans = {"FT9","FT10","TP9","TP10"} % Noisy_Channels
    cfg.f1_grid = 0.1:0.1:0.5
    cfg.f2_grid = 1:1:100
    cfg.orders  = 2:10
    cfg.save_dir (1,1) string = "results/train_results"
    cfg.Fs (1,1) double = 500
    cfg.target_group (1,1) string = "group1"
    cfg.kfold = "loocv"
    cfg.is_norm_proj (1,1) double = 0
    cfg.is_preprocessed (1,1) logical = false
    cfg.detect_additional_noisy (1,1) logical = true
end

% Create save directory if it doesn't exist
if ~exist(cfg.save_dir,'dir')
    mkdir(cfg.save_dir);
end

% ---- Load training data initially to check preprocessing status ----
fprintf('\n Loading training data from: %s\n', cfg.data_train);
[DataTrain, ChannelLoc, SubjectIDs, GroupNames] = utils.load_data(cfg.data_train, cfg.channel_list, {});

% Auto-detect if preprocessing was already done if not explicitly set
if ~cfg.is_preprocessed
    cfg.is_preprocessed = utils.detect_preprocessing_status(DataTrain);
end

fprintf('\n Training mode: %s\n', cfg.mode);
fprintf('Preprocessing status: %s\n', string(cfg.is_preprocessed));

% Handle noisy channels based on preprocessing status
if cfg.is_preprocessed
    fprintf('Preprocessing already done - skipping noisy channel exclusion\n');
    effective_exclude_chans = {}; % No channels to exclude
else
    fprintf('Preprocessing not done - will exclude predefined noisy channels: %s\n', ...
        strjoin(cfg.exclude_chans, ', '));
    effective_exclude_chans = cfg.exclude_chans;
    
    % Additional noisy channel detection if requested
    if cfg.detect_additional_noisy
        fprintf('Detecting additional noisy channels...\n');
        additional_noisy = utils.detect_additional_noisy_channels(DataTrain);
        if ~isempty(additional_noisy)
            effective_exclude_chans = union(effective_exclude_chans, additional_noisy);
            fprintf('Additional noisy channels detected: %s\n', strjoin(additional_noisy, ', '));
        else
            fprintf('No additional noisy channels detected\n');
        end
    end
end

% Reload data with appropriate channel exclusion
if ~isempty(effective_exclude_chans)
    fprintf('Excluding channels: %s\n', strjoin(effective_exclude_chans, ', '));
    [DataTrain, ChannelLoc, SubjectIDs, GroupNames] = utils.load_data(cfg.data_train, cfg.channel_list, effective_exclude_chans);
end

channels = keys(DataTrain);
Nch = numel(channels);

fprintf('Channels available for analysis: %d\n', Nch);
if Nch == 0
    error('No channels available for analysis after exclusion');
end

% Initialize results structure
BestParamsAll = struct('channel',{}, 'cutoff',{}, 'order',{}, 'dim',{}, ...
    'metric',{}, 'mode',{}, 'polarity',{}, 'GroupNames',{}, 'Fs',{}, ...
    'P0_train',{}, 'm0_train',{}, 'P1_train',{}, 'm1_train',{}, ...
    'is_preprocessed',{});

% Pre-read labels if needed for correlation
labelsMap = containers.Map();
if cfg.mode == "correlation"
    if isempty(cfg.labels_file) || ~isfile(cfg.labels_file)
        error('For correlation mode, labels_file must be provided and exist');
    end
    fprintf('Loading clinical labels from: %s\n', cfg.labels_file);
    T = utils.read_labels_table(cfg.labels_file);
    if isempty(T)
        error('Could not load labels from file');
    end
    ids = string(T.ID);
    y   = double(T.Target);
    for i = 1:numel(ids)
        labelsMap(ids(i)) = y(i);
    end
    fprintf('Loaded %d subject labels\n', numel(ids));
end

% Display grid search information
fprintf('\n Starting grid search:\n');
fprintf('   Frequency grid: f1=%.1f:%.1f:%.1f, f2=%d:%d:%d\n', ...
    cfg.f1_grid(1), cfg.f1_grid(2)-cfg.f1_grid(1), cfg.f1_grid(end), ...
    cfg.f2_grid(1), cfg.f2_grid(2)-cfg.f2_grid(1), cfg.f2_grid(end));
fprintf('   AR orders: %d:%d:%d\n', cfg.orders(1), cfg.orders(2)-cfg.orders(1), cfg.orders(end));
fprintf('   Using %s cross-validation\n', cfg.kfold);

% Start parallel processing if available
if isempty(gcp('nocreate'))
    try
        parpool;
        fprintf('Starting parallel pool...\n');
    catch
        fprintf('Parallel pool not available, running sequentially\n');
    end
end

% Main training loop per channel
fprintf('\n🎯 Training per channel:\n');
parfor chIdx = 1:Nch
    ch = channels{chIdx};
    chData = DataTrain(ch);
    
    fprintf('   Processing channel: %s\n', ch);
    
    % Assemble AllData and Classes
    if isfield(chData, 'group2') && ~isempty(chData.group2)
        % Two-group data
        AllData = [chData.group1, chData.group2];
        Classes = [ones(numel(chData.group1), 1); zeros(numel(chData.group2), 1)];
        n1 = numel(chData.group1);
        n2 = numel(chData.group2);
        fprintf('     Two groups: %d (Group1) + %d (Group2) = %d total subjects\n', n1, n2, n1+n2);
    else
        % One-group data
        AllData = chData.group1;
        Classes = ones(numel(AllData), 1);
        n1 = numel(AllData);
        n2 = 0;
        fprintf('     Single group: %d subjects\n', n1);
    end

    bestMetric = -inf;
    bestRec = struct();
    total_combinations = numel(cfg.f1_grid) * numel(cfg.f2_grid) * numel(cfg.orders);
    processed_combinations = 0;

    for f1 = cfg.f1_grid
        for f2 = cfg.f2_grid
            % Skip if frequency band is too narrow
            if (f2 - f1) < 4
                continue;
            end
            
            Filt = [0 f1 0; f2 inf 0];
            
            % Filter all subjects for current frequency band
            Xf = utils.filter_data(AllData, Filt, cfg.Fs);
            
            for ord = cfg.orders
                processed_combinations = processed_combinations + 1;
                
                % Compute LPC coefficients for all subjects
                LPC_all = utils.compute_yw(Xf, ord);
                dims = 1:(ord-1);
                
                if cfg.mode == "classification"
                    % ---- CLASSIFICATION TRAINING ----
                    for d = dims
                        acc = local_loocv_accuracy(LPC_all, Classes, d, cfg.is_norm_proj);
                        
                        if acc > bestMetric
                            bestMetric = acc;
                            
                            % Store best parameters
                            bestRec.channel = ch;
                            bestRec.cutoff = [f1 f2];
                            bestRec.order = ord;
                            bestRec.dim = d;
                            bestRec.metric = acc;
                            bestRec.mode = 'classification';
                            bestRec.polarity = 0;
                            
                            % SAVE HYPERPLANES FROM TRAINING DATA
                            if n2 > 0  % Two groups
                                [bestRec.P0_train, bestRec.m0_train] = utils.build_hyperplanes(LPC_all(n1+1:end,:), d);
                                [bestRec.P1_train, bestRec.m1_train] = utils.build_hyperplanes(LPC_all(1:n1,:), d);
                            else  % One group (use same data for both)
                                [bestRec.P0_train, bestRec.m0_train] = utils.build_hyperplanes(LPC_all, d);
                                [bestRec.P1_train, bestRec.m1_train] = utils.build_hyperplanes(LPC_all, d);
                            end
                        end
                    end
                    
                else
                    % ---- CORRELATION TRAINING ----
                    % Determine target group indices
                    if cfg.target_group == "group1"
                        tgt_ids = SubjectIDs.(ch).group1;
                        tgt_idx = 1:numel(tgt_ids);
                    else
                        tgt_ids = SubjectIDs.(ch).group2;
                        tgt_idx = (n1+1):(n1+numel(chData.group2));
                    end
                    
                    % Get target values
                    y = utils.fetch_targets(tgt_ids, labelsMap);
                    if all(isnan(y)) || numel(y) < 3
                        continue;
                    end
                    
                    for d = dims
                        % Build reference hyperplane from non-target group
                        if isfield(chData, 'group2') && ~isempty(chData.group2)
                            if cfg.target_group == "group1"
                                [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(n1+1:end,:), d);
                                [P_tar_full, m_tar_full] = utils.build_hyperplanes(LPC_all(1:n1,:), d);
                            else
                                [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(1:n1,:), d);
                                [P_tar_full, m_tar_full] = utils.build_hyperplanes(LPC_all(n1+1:end,:), d);
                            end
                        else
                            % Single group - use same data for reference and target
                            [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(tgt_idx,:), d);
                            [P_tar_full, m_tar_full] = utils.build_hyperplanes(LPC_all(tgt_idx,:), d);
                        end

                        % LOOCV within target group
                        scores = zeros(numel(tgt_idx), 1);
                        for i = 1:numel(tgt_idx)
                            ti = tgt_idx(i);
                            if cfg.target_group == "group1"
                                tr_idx = setdiff(1:n1, ti);
                                [P_tar, m_tar] = utils.build_hyperplanes(LPC_all(tr_idx,:), d);
                            else
                                base = n1;
                                loc = ti - base;
                                tr_idx = base + setdiff(1:numel(chData.group2), loc);
                                [P_tar, m_tar] = utils.build_hyperplanes(LPC_all(tr_idx,:), d);
                            end
                            scores(i) = utils.compute_leapd_scores(LPC_all(ti,:), P_ref, m_ref, P_tar, m_tar, d, cfg.is_norm_proj);
                        end

                        % Find optimal polarity
                        [rho, p, pol, scores_aligned] = utils.pick_polarity_and_rho(scores, y);
                        
                        if abs(rho) > abs(bestMetric)
                            bestMetric = rho;
                            
                            % Store best parameters
                            bestRec.channel = ch;
                            bestRec.cutoff = [f1 f2];
                            bestRec.order = ord;
                            bestRec.dim = d;
                            bestRec.metric = rho;
                            bestRec.mode = 'correlation';
                            bestRec.polarity = pol;
                            
                            % SAVE HYPERPLANES FROM TRAINING DATA
                            bestRec.P0_train = P_ref;
                            bestRec.m0_train = m_ref;
                            bestRec.P1_train = P_tar_full;  % Use full target group for testing
                            bestRec.m1_train = m_tar_full;
                        end
                    end
                end
            end
        end
    end

    % Store best results for this channel
    if ~isempty(fieldnames(bestRec))
        bestRec.GroupNames = GroupNames;
        bestRec.Fs = cfg.Fs;
        bestRec.is_preprocessed = cfg.is_preprocessed;
        BestParamsAll(chIdx) = bestRec;
        
        if cfg.mode == "classification"
            fprintf('Best: f=[%.1f-%.1f]Hz, order=%d, dim=%d, ACC=%.2f%%\n', ...
                bestRec.cutoff(1), bestRec.cutoff(2), bestRec.order, bestRec.dim, bestRec.metric);
        else
            fprintf('Best: f=[%.1f-%.1f]Hz, order=%d, dim=%d, ρ=%.3f (p=%.3f)\n', ...
                bestRec.cutoff(1), bestRec.cutoff(2), bestRec.order, bestRec.dim, bestRec.metric, p);
        end
    else
        fprintf('No valid parameters found for channel %s\n', ch);
    end
end

% Remove empty entries
empty_mask = arrayfun(@(x) isempty(x.channel), BestParamsAll);
BestParamsAll = BestParamsAll(~empty_mask);

% Save results
metadata = struct('timestamp', datestr(now), 'cfg', cfg, 'is_preprocessed', cfg.is_preprocessed);
save(fullfile(cfg.save_dir, 'BestParamsAll.mat'), 'BestParamsAll', 'metadata');

results.BestParamsAll = BestParamsAll;
results.metadata = metadata;

fprintf('\n Training completed successfully!\n');
fprintf('Results saved to: %s\n', fullfile(cfg.save_dir, 'BestParamsAll.mat'));
fprintf('Total channels with valid models: %d/%d\n', numel(BestParamsAll), Nch);
fprintf('Preprocessing status: %s\n', string(cfg.is_preprocessed));

end

% Local function for LOOCV accuracy calculation
function acc = local_loocv_accuracy(LPC_all, Classes, d, is_norm)
% Leave-one-out cross-validation accuracy for classification

    if numel(unique(Classes)) < 2
        acc = NaN;
        return;
    end
    
    n = numel(Classes);
    preds = false(n, 1);
    
    for i = 1:n
        tr = setdiff(1:n, i);
        c1 = find(Classes(tr) == 1);
        c0 = find(Classes(tr) == 0);
        
        [P1, m1] = utils.build_hyperplanes(LPC_all(tr(c1), :), d);
        [P0, m0] = utils.build_hyperplanes(LPC_all(tr(c0), :), d);
        
        s = utils.compute_leapd_scores(LPC_all(i, :), P0, m0, P1, m1, d, is_norm);
        preds(i) = s >= 0.5;  % CORRECT: >=0.5 = Group 1, <0.5 = Group 2
    end
    
    acc = mean(double(preds) == Classes) * 100;
end
