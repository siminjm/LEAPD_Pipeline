function results = main_train(cfg)
% MAIN_TRAIN - Unified training for classification and correlation.

% cfg fields:
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

% ---- Load training data ----
[DataTrain, ChannelLoc, SubjectIDs, GroupNames] = utils.load_data(cfg.data_train, cfg.channel_list, cfg.exclude_chans);
channels = keys(DataTrain);
Nch = numel(channels);

fprintf('\n Training mode: %s | channels: %d\n', cfg.mode, Nch);

BestParamsAll = struct('channel',{},'cutoff',{},'order',{},'dim',{},...
    'metric',{},'mode',{},'polarity',{},'GroupNames',{},'Fs',{},...
    'P0_train',{},'m0_train',{},'P1_train',{},'m1_train',{}); % ADDED: Save hyperplanes

% Pre-read labels if needed
labelsMap = containers.Map();
if cfg.mode == "correlation"
    T = utils.read_labels_table(cfg.labels_file);
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
        n1 = numel(chData.group1);
        n2 = numel(chData.group2);
    else
        % one-group
        AllData = chData.group1;
        Classes = ones(numel(AllData),1);
        n1 = numel(AllData);
        n2 = 0;
    end

    bestMetric = -inf; 
    bestRec = struct();

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
                            bestMetric = acc; 
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
                    % ---- Correlation training
                    if cfg.target_group == "group1"
                        tgt_ids = SubjectIDs.(ch).group1; 
                        tgt_idx = 1:numel(tgt_ids);
                    else
                        tgt_ids = SubjectIDs.(ch).group2;
                        tgt_idx = (n1+1):(n1+numel(chData.group2));
                    end
                    y = utils.fetch_targets(tgt_ids, labelsMap);
                    if all(isnan(y)) || numel(y) < 3, continue; end

                    dims = 1:(ord-1);
                    for d = dims
                        % Build reference from non-target group
                        if isfield(chData,'group2') && ~isempty(chData.group2)
                            if cfg.target_group == "group1"
                                [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(n1+1:end,:), d);
                                [P_tar_full, m_tar_full] = utils.build_hyperplanes(LPC_all(1:n1,:), d);
                            else
                                [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(1:n1,:), d);
                                [P_tar_full, m_tar_full] = utils.build_hyperplanes(LPC_all(n1+1:end,:), d);
                            end
                        else
                            % pseudo-reference from target itself
                            [P_ref, m_ref] = utils.build_hyperplanes(LPC_all(tgt_idx,:), d);
                            [P_tar_full, m_tar_full] = utils.build_hyperplanes(LPC_all(tgt_idx,:), d);
                        end

                        % LOOCV within target group
                        scores = zeros(numel(tgt_idx),1);
                        for i = 1:numel(tgt_idx)
                            ti = tgt_idx(i);
                            if cfg.target_group == "group1"
                                tr_idx = setdiff(1:n1, ti);
                                [P_tar, m_tar] = utils.build_hyperplanes(LPC_all(tr_idx,:), d);
                            else
                                base = n1;
                                loc  = ti - base;
                                tr_idx = base + setdiff(1:numel(chData.group2), loc);
                                [P_tar, m_tar] = utils.build_hyperplanes(LPC_all(tr_idx,:), d);
                            end
                            scores(i) = utils.compute_leapd_scores(LPC_all(ti,:), P_ref, m_ref, P_tar, m_tar, d);
                        end

                        [rho, p, pol, scores_aligned] = utils.pick_polarity_and_rho(scores, y);
                        if abs(rho) > abs(bestMetric)
                            bestMetric = rho; 
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

    if ~isempty(fieldnames(bestRec))
        bestRec.GroupNames = GroupNames; 
        bestRec.Fs = cfg.Fs;
        BestParamsAll(chIdx) = bestRec;
    end
end

metadata = struct('timestamp',datestr(now), 'cfg',cfg);
save(fullfile(cfg.save_dir,'BestParamsAll.mat'), 'BestParamsAll','metadata');
results.BestParamsAll = BestParamsAll; 
results.metadata = metadata;
fprintf('\n Saved training results to %s\n', fullfile(cfg.save_dir,'BestParamsAll.mat'));
end

function acc = local_loocv_accuracy(LPC_all, Classes, d, is_norm)
% Leave-one-out across all subjects
if numel(unique(Classes))<2, acc = NaN; return; end
n = numel(Classes); 
preds = false(n,1);
for i=1:n
    tr = setdiff(1:n,i);
    c1 = find(Classes(tr)==1);
    c0 = find(Classes(tr)==0);
    [P1,m1] = utils.build_hyperplanes(LPC_all(tr(c1),:), d);
    [P0,m0] = utils.build_hyperplanes(LPC_all(tr(c0),:), d);
    s = utils.compute_leapd_scores(LPC_all(i,:), P0,m0, P1,m1, d, is_norm);
    preds(i) = s >= 0.5;  % CORRECT: >=0.5 = Group 1
end
acc = mean(double(preds)==Classes)*100;
end
