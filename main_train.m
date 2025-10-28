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

% Outputs (also saved as MAT under save_dir/BestParamsAll.mat):
%   results.BestParamsAll  : struct array with per-channel best params
%   results.metadata       : cfg + timestamp

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
end

if ~exist(cfg.save_dir,'dir'); mkdir(cfg.save_dir); end

% ---- Load training data (generic 1-or-2 group) ----
[DataTrain, ChannelLoc, SubjectIDs, GroupNames] = utils.load_data(cfg.data_train, cfg.channel_list, cfg.exclude_chans);
% DataTrain is a containers.Map from channel-> struct('group1',{cell subjects}, 'group2',{cell subjects or []})

channels = keys(DataTrain);
Nch = numel(channels);

fprintf('\n🚀 Training mode: %s | channels: %d\n', cfg.mode, Nch);

BestParamsAll = struct('channel',{},'cutoff',{},'order',{},'dim',{},...
    'metric',{},'mode',{},'polarity',{},'GroupNames',{},'Fs',{});

% Pre-read labels if needed
labelsMap = containers.Map();
if cfg.mode == "correlation"
    T = utils.read_labels_table(cfg.labels_file);
    % Expect columns: ID, Target
    ids = string(T.ID);
    y   = double(T.Target);
    for i=1:numel(ids), labelsMap(ids(i)) = y(i); end
end

parfor chIdx = 1:Nch
    ch = channels{chIdx};
    chData = DataTrain(ch);

    % Assemble AllData for filtering shape consistency
    if isfield(chData,'group2') && ~isempty(chData.group2)
        % two-group
        AllData = [chData.group1, chData.group2];
        Classes = [ones(numel(chData.group1),1); zeros(numel(chData.group2),1)];
    else
        % one-group
        AllData = chData.group1;
        Classes = ones(numel(AllData),1);
    end

    bestMetric = -inf; bestRec = struct();

    for f1 = cfg.f1_grid
        for f2 = cfg.f2_grid
            if (f2 - f1) < 4, continue; end
            Filt = [0 f1 0; f2 inf 0];
            % Filter all subjects once for speed
            Xf = utils.filter_data(AllData, Filt, cfg.Fs);
            % Precompute LPC per order
            for ord = cfg.orders
                LPC_all = utils.compute_yw(Xf, ord);

                if cfg.mode == "classification"
                    % ---- Classification: LOOCV across both groups
                    dims = 1:(ord-1);
                    for d = dims
                        acc = local_loocv_accuracy(LPC_all, Classes, d, cfg.is_norm_proj);
                        if acc > bestMetric
                            bestMetric = acc; bestRec.channel = ch; bestRec.cutoff = [f1 f2];
                            bestRec.order = ord; bestRec.dim = d; bestRec.metric = acc;
                            bestRec.mode = 'classification'; bestRec.polarity = 0; %#ok<PFBNS>
                        end
                    end
                else
                    % ---- Correlation training
                    % Determine target vector over the target group
                    if cfg.target_group == "group1"
                        tgt_ids = SubjectIDs.(ch).group1; tgt_idx = 1:numel(tgt_ids);
                    else
                        tgt_ids = SubjectIDs.(ch).group2;
                        tgt_idx = (numel(chData.group1)+1) : (numel(chData.group1)+numel(chData.group2));
                    end
                    y = utils.fetch_targets(tgt_ids, labelsMap);
                    if all(isnan(y)) || numel(y) < 3, continue; end

                    dims = 1:(ord-1);
                    for d = dims
                        % Build reference once: non-target group if exists; otherwise, pseudo-reference
                        if isfield(chData,'group2') && ~isempty(chData.group2)
                            if cfg.target_group == "group1"
                                [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(numel(chData.group1)+1:end,:), d);
                            else
                                [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(1:numel(chData.group1),:), d);
                            end
                        else
                            % pseudo-reference from target itself (mean-only plane)
                            [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(tgt_idx,:), d);
                        end

                        % LOOCV within target group: compare distances (ref vs. target plane)
                        scores = zeros(numel(tgt_idx),1);
                        for i = 1:numel(tgt_idx)
                            ti = tgt_idx(i);
                            if cfg.target_group == "group1"
                                tr_idx = setdiff(1:numel(chData.group1), ti);
                                [P_tar, m_tar] = utils.build_hyperplanes(LPC_all(tr_idx,:), d);
                            else
                                base = numel(chData.group1);
                                loc  = ti - base; % index within group2
                                tr_idx = base + setdiff(1:numel(chData.group2), loc);
                                [P_tar, m_tar] = utils.build_hyperplanes(LPC_all(tr_idx,:), d);
                            end
                            scores(i) = utils.compute_leapd_scores(LPC_all(ti,:), P_ref, m_ref, P_tar, m_tar, 1); % one subject vector
                        end

                        [rho, p, pol, ~] = utils.pick_polarity_and_rho(scores, y);
                        if abs(rho) > abs(bestMetric)
                            bestMetric = rho; bestRec.channel = ch; bestRec.cutoff = [f1 f2];
                            bestRec.order = ord; bestRec.dim = d; bestRec.metric = rho;
                            bestRec.mode = 'correlation'; bestRec.polarity = pol; %#ok<PFBNS>
                        end
                    end
                end
            end
        end
    end

    if ~isempty(fieldnames(bestRec))
        bestRec.GroupNames = GroupNames; bestRec.Fs = cfg.Fs;
        BestParamsAll(chIdx) = bestRec; %#ok<PFOUS>
    end
end

metadata = struct('timestamp',datestr(now), 'cfg',cfg);
save(fullfile(cfg.save_dir,'BestParamsAll.mat'), 'BestParamsAll','metadata');
results.BestParamsAll = BestParamsAll; results.metadata = metadata;
fprintf('\n Saved training results to %s\n', fullfile(cfg.save_dir,'BestParamsAll.mat'));
end

function acc = local_loocv_accuracy(LPC_all, Classes, d, is_norm)
% Leave-one-out across all subjects (two-class only; if single-class, return NaN)
if numel(unique(Classes))<2, acc = NaN; return; end
n = numel(Classes); preds = false(n,1);
for i=1:n
    tr = setdiff(1:n,i);
    yte = Classes(i); %#ok<NASGU>
    c1 = find(Classes(tr)==1);
    c0 = find(Classes(tr)==0);
    [P1,m1] = utils.build_hyperplanes(LPC_all(tr(c1),:), d);
    [P0,m0] = utils.build_hyperplanes(LPC_all(tr(c0),:), d);
    s = utils.compute_leapd_scores(LPC_all(i,:), P0,m0, P1,m1, 1, is_norm);
    preds(i) = s >= 0.5;
end
acc = mean(double(preds)==Classes)*100;
end
