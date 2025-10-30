function results = main_test(cfg)
% MAIN_TEST - Unified test for classification and correlation.

arguments
    cfg.mode (1,1) string {mustBeMember(cfg.mode,["classification","correlation"])}
    cfg.data_test (1,1) string
    cfg.trained_model (1,1) string
    cfg.labels_file (1,1) string = ""
    cfg.combo_sizes = 1:10
    cfg.max_full_combos (1,1) double = 5
    cfg.Fs (1,1) double = 500
    cfg.is_norm_proj (1,1) double = 0
    cfg.save_dir (1,1) string = "results/test_results"
end

if ~exist(cfg.save_dir,'dir'); mkdir(cfg.save_dir); end

S = load(cfg.trained_model,'BestParamsAll'); 
BestParamsAll = S.BestParamsAll;
[DataTest, ~, SubjectIDs, ~] = utils.load_data(cfg.data_test, [], {});

% Keep only channels that exist in both train results and test data
trained_chs = string({BestParamsAll.channel});
channels_all = keys(DataTest);
avail = ismember(trained_chs, string(channels_all));
BestParamsAll = BestParamsAll(avail);
channels = string({BestParamsAll.channel});
Nch = numel(channels);

fprintf('\n Testing mode: %s | channels usable: %d\n', cfg.mode, Nch);

% ---- Precompute single-channel scores on test subjects ----
totalSubs = utils.count_subjects(DataTest);
ScoresMatrix = NaN(totalSubs, Nch);
TrueLabels = [];
TargetVec   = [];

% ---- Correlation labels (if applicable) ----
labelsMap = containers.Map();
if cfg.mode == "correlation" && strlength(cfg.labels_file) > 0
    T = utils.read_labels_table(cfg.labels_file);
    if ~isempty(T)
        ids = string(T.ID);
        y = double(T.Target);
        for i = 1:numel(ids)
            labelsMap(ids(i)) = y(i);
        end
    end
end

rowOffset = 0;
chan_index_map = strings(Nch,1);

for chIdx = 1:Nch
    ch = channels(chIdx);
    P = BestParamsAll(chIdx);
    chData = DataTest(char(ch));

    % --- Build test subjects list and labels ---
    if isfield(chData,'group2') && ~isempty(chData.group2)
        AllData = [chData.group1, chData.group2];
        TL = [ones(numel(chData.group1),1); zeros(numel(chData.group2),1)];
        
        if cfg.mode == "correlation"
            grp_ids = [SubjectIDs.(char(ch)).group1; SubjectIDs.(char(ch)).group2];
            TargetVec = utils.fetch_targets(grp_ids, labelsMap);
        end
    else
        AllData = chData.group1;
        TL = ones(numel(AllData),1);
        
        if cfg.mode == "correlation"
            grp_ids = SubjectIDs.(char(ch)).group1;
            TargetVec = utils.fetch_targets(grp_ids, labelsMap);
        end
    end

    if isempty(TrueLabels)
        TrueLabels = TL;
    end

    % --- APPLY TRAINED PARAMETERS TO TEST DATA ---
    % 1. Use trained filter cutoffs
    Filt = [0 P.cutoff(1) 0; P.cutoff(2) inf 0];
    Xf = utils.filter_data(AllData, Filt, cfg.Fs);
    
    % 2. Use trained AR order
    LPC_test = utils.compute_yw(Xf, P.order);
    
    % 3. Use trained hyperplanes from TRAINING data
    s = utils.compute_leapd_scores(LPC_test, P.P0_train, P.m0_train, P.P1_train, P.m1_train, P.dim, cfg.is_norm_proj);

    % 4. Apply trained polarity for correlation
    if cfg.mode=="correlation" && isfield(P,'polarity') && P.polarity==-1
        s = 1 - s;
    end

    ScoresMatrix(1:numel(s), chIdx) = s(:);
    chan_index_map(chIdx) = ch;
end

% ---- Evaluate combinations (same as before) ----
results = struct();
results.single = struct();
results.combos = struct();

% Single-channel metrics
for j=1:Nch
    if cfg.mode=="classification"
        metrics = utils.evaluate_classification(ScoresMatrix(:,j), TrueLabels);
    else
        metrics = utils.evaluate_correlation(ScoresMatrix(:,j), TargetVec);
    end
    results.single(j).channel = chan_index_map(j);
    results.single(j).metrics = metrics;
end

best_prev = struct();
for k = cfg.combo_sizes
    if k <= cfg.max_full_combos
        combos = nchoosek(1:Nch, k);
    else
        combos = utils.generate_combinations(best_prev, Nch, k);
    end

    best_k = struct('indices',[],'channels',[],'metrics',[],'score',-inf);
    for i=1:size(combos,1)
        idx = combos(i,:);
        s = utils.combine_scores(ScoresMatrix(:,idx));
        if cfg.mode=="classification"
            metrics = utils.evaluate_classification(s, TrueLabels);
            score = metrics.ACC;
        else
            metrics = utils.evaluate_correlation(s, TargetVec);
            score = abs(metrics.Rho);
        end
        if score > best_k.score
            best_k.indices = idx;
            best_k.channels = cellstr(chan_index_map(idx));
            best_k.metrics = metrics;
            best_k.score = score;
        end
    end
    results.combos(k).best = best_k;
    if k <= cfg.max_full_combos
        best_prev(k).indices = best_k.indices;
    end
    fprintf('k=%d best: %s | ', k, strjoin(best_k.channels,', '));
    if cfg.mode=="classification"
        fprintf('ACC=%.2f%%\n', best_k.metrics.ACC);
    else
        fprintf('Rho=%.3f (p=%.3g)\n', best_k.metrics.Rho, best_k.metrics.PValue);
    end
end

save(fullfile(cfg.save_dir,'test_results.mat'),'results','cfg');
fprintf('\n Saved test results to %s\n', fullfile(cfg.save_dir,'test_results.mat'));
end
