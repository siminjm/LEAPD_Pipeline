function results = main_test(cfg)
% MAIN_TEST - Unified test for classification and correlation.

% cfg fields:
%   .mode             : "classification" | "correlation"
%   .data_test        : path to test EEG .mat
%   .trained_model    : path to results/train_results/BestParamsAll.mat
%   .labels_file      : target file (for correlation testing)
%   .combo_sizes      : vector (e.g., 1:10)
%   .max_full_combos  : integer (default 5) – enumerate up to this, then use fixed-best
%   .Fs               : sampling rate (default 500)
%   .is_norm_proj     : 0 or 1
%   .save_dir         : 'results/test_results'

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

S = load(cfg.trained_model,'BestParamsAll'); BestParamsAll = S.BestParamsAll;
[DataTest, ~, SubjectIDs, ~] = utils.load_data(cfg.data_test, [], {});
channels_all = keys(DataTest);

% Keep only channels that exist in both train results and test data
trained_chs = string({BestParamsAll.channel});
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
labelsMap = containers.Map();   % initialize empty map
T = table();                    % predefine to avoid scope errors

if cfg.mode == "correlation" && strlength(cfg.labels_file) > 0
    fprintf('\n Loading correlation targets from %s...\n', cfg.labels_file);

    T = utils.read_labels_table(cfg.labels_file);
    if isempty(T)
        warning('No label table loaded; correlation metrics will be skipped.');
    else
        ids = string(T.ID);
        y = double(T.Target);
        for i = 1:numel(ids)
            labelsMap(ids(i)) = y(i);
        end
    end
end


rowOffset = 0; % fill by first channel size (assumes consistent subject ordering)
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
        TargetVec = utils.fetch_targets(grp_ids, labelsMap); % allow NaNs
    end

else
    AllData = chData.group1;
    TL = ones(numel(AllData),1);

    if cfg.mode == "correlation"
        grp_ids = SubjectIDs.(char(ch)).group1;
        TargetVec = utils.fetch_targets(grp_ids, labelsMap);
    end
end

% --- Safety check for missing targets ---
if cfg.mode == "correlation" && sum(isnan(TargetVec)) > 0
    missing = grp_ids(isnan(TargetVec));
    warn_subset = missing(1:min(5,numel(missing))); % show up to 5 IDs
    warning('%d of %d subject IDs in test set had no matching entry in %s.\nMissing examples: %s', ...
        sum(isnan(TargetVec)), numel(TargetVec), cfg.labels_file, strjoin(string(warn_subset), ', '));
end

if isempty(TrueLabels)
    TrueLabels = TL;
end


    % filter + LPC using trained hyperparams
    Filt = [0 P.cutoff(1) 0; P.cutoff(2) inf 0];
    Xf = utils.filter_data(AllData, Filt, cfg.Fs);
    LPC_all = utils.compute_yw(Xf, P.order);

    % Build planes
    if isfield(chData,'group2') && ~isempty(chData.group2)
        n1 = numel(chData.group1);
        [P0,m0] = utils.build_hyperplanes(LPC_all(n1+1:end,:), P.dim); % class 0 ref
        [P1,m1] = utils.build_hyperplanes(LPC_all(1:n1,:), P.dim);     % class 1 tar
    else
        [P0,m0] = utils.build_hyperplanes(LPC_all, P.dim);
        [P1,m1] = utils.build_hyperplanes(LPC_all, P.dim);
    end

    s = utils.compute_leapd_scores(LPC_all, P0,m0, P1,m1, [], cfg.is_norm_proj);

    % Polarity alignment for correlation
    if cfg.mode=="correlation" && isfield(P,'polarity') && P.polarity==-1
        s = 1 - s;
    end

    ScoresMatrix(1:numel(s),chIdx) = s(:);
    chan_index_map(chIdx) = ch;
end

% ---- Evaluate combinations ----
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

best_prev = struct(); % hold best combos for k<=5 to seed >5

for k = cfg.combo_sizes
    if k <= cfg.max_full_combos
        combos = nchoosek(1:Nch, k);
    else
        combos = utils.generate_combinations(best_prev, Nch, k); % progressive fixed-best
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
        best_prev(k).indices = best_k.indices; %#ok<AGROW>
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
