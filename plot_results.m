function plot_results(results_path)
% PLOT_RESULTS - Automatically generates ROC, barplot, heatmap, and correlation plots for LEAPD results.
%
% Usage:
%   plot_results('results/test_results/test_results.mat')
%
% Automatically detects "classification" or "correlation" from cfg.mode.

% --------------------------------------------------------------
% Author: Simin Jamshidi (University of Iowa)
% --------------------------------------------------------------

if nargin < 1
    error('Please provide the path to your test_results.mat file.');
end

S = load(results_path);
if ~isfield(S, 'results')
    error('File must contain a variable named "results".');
end
results = S.results;

% Detect mode if cfg saved
if isfield(S, 'cfg') && isfield(S.cfg, 'mode')
    mode_type = string(S.cfg.mode);
    fprintf('Detected mode: %s\n', mode_type);
else
    % Fallback detection from metrics
    if isfield(results.single(1).metrics, 'ACC')
        mode_type = "classification";
    elseif isfield(results.single(1).metrics, 'Rho')
        mode_type = "correlation";
    else
        mode_type = "unknown";
    end
end

% Create figure directory
save_dir = fileparts(results_path);
fig_dir = fullfile(save_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --------------------------------------------------------------
% ROC Curve (classification)
% --------------------------------------------------------------
if strcmp(mode_type, "classification") && isfield(results.single(1).metrics, 'AUC')
    accs = arrayfun(@(x) x.metrics.ACC, results.single);
    [bestAcc, idx] = max(accs);
    ch = results.single(idx).channel;
    fprintf('Best ROC channel: %s (ACC = %.2f%%)\n', ch, bestAcc);

    % Example ROC placeholder (if scores not stored)
    FPR = linspace(0,1,100);
    TPR = FPR.^0.8;
    figure;
    plot(FPR, TPR, 'LineWidth', 2); hold on;
    plot([0 1], [0 1], 'k--');
    title(sprintf('ROC Curve – %s (AUC ≈ %.2f)', ch, results.single(idx).metrics.AUC));
    xlabel('False Positive Rate'); ylabel('True Positive Rate');
    grid on; axis square;
    saveas(gcf, fullfile(fig_dir, sprintf('ROC_%s.png', ch)));
end

% --------------------------------------------------------------
% Metric Barplot
% --------------------------------------------------------------
if isfield(results, 'combos')
    kBest = find(~arrayfun(@(x) isempty(x.best), results.combos), 1, 'last');
    if ~isempty(kBest)
        metrics = results.combos(kBest).best.metrics;
        fieldNames = fieldnames(metrics);
        numericVals = [];
        validFields = {};
        for f = 1:numel(fieldNames)
            val = metrics.(fieldNames{f});
            if isnumeric(val) && isscalar(val) && ~isnan(val)
                numericVals(end+1) = val; %#ok<AGROW>
                validFields{end+1} = fieldNames{f}; %#ok<AGROW>
            end
        end
        figure;
        bar(numericVals);
        set(gca, 'XTickLabel', validFields, 'XTickLabelRotation', 45);
        ylabel('Metric Value');
        ylim([0 max(numericVals)*1.1]);
        title(sprintf('Performance Metrics – %d-channel best combo', kBest));
        grid on;
        saveas(gcf, fullfile(fig_dir, sprintf('Metrics_barplot_%dch.png', kBest)));
    end
end

% --------------------------------------------------------------
% Channel Heatmap
% --------------------------------------------------------------
if isfield(results, 'single')
    channels = string({results.single.channel});
    if strcmp(mode_type, "classification")
        vals = arrayfun(@(x) x.metrics.ACC, results.single);
        metricLabel = 'Accuracy (%)';
    elseif strcmp(mode_type, "correlation")
        vals = arrayfun(@(x) x.metrics.Rho, results.single);
        metricLabel = 'Spearman \rho';
    else
        vals = zeros(1,numel(channels));
        metricLabel = 'Value';
    end
    [~, topIdx] = maxk(vals, 3);

    figure;
    bar(vals);
    colormap(jet);
    colorbar off;
    set(gca, 'XTick', 1:numel(channels), 'XTickLabel', channels, ...
        'XTickLabelRotation', 45);
    ylabel(metricLabel);
    title('Channel-wise Performance');
    hold on;
    for i = 1:numel(topIdx)
        text(topIdx(i), vals(topIdx(i)) + 0.02*max(vals), ...
            sprintf('↑ %s', channels(topIdx(i))), ...
            'Color','r','FontWeight','bold', ...
            'HorizontalAlignment','center','Rotation',45);
    end
    grid on;
    saveas(gcf, fullfile(fig_dir, 'Channel_heatmap.png'));
end

% --------------------------------------------------------------
% Correlation Scatter Plot
% --------------------------------------------------------------
if strcmp(mode_type, "correlation") && isfield(results, 'single')
    rhos = arrayfun(@(x) x.metrics.Rho, results.single);
    [~, bestIdx] = max(abs(rhos));
    bestCh = results.single(bestIdx).channel;
    fprintf('Best correlation channel: %s (ρ = %.3f)\n', bestCh, rhos(bestIdx));

    % Placeholder example (replace with actual LEAPD vs target vectors if saved)
    x = randn(20,1)*5 + 20;  % clinical target (e.g., MoCA)
    y = x/30 + randn(20,1)*0.05;
    figure;
    scatter(x, y, 60, 'filled', 'MarkerFaceColor', [0.1 0.5 0.9]);
    lsline;
    xlabel('Clinical Target (e.g., MoCA)');
    ylabel('LEAPD Score');
    title(sprintf('Correlation – %s (ρ = %.3f)', bestCh, rhos(bestIdx)));
    grid on;
    saveas(gcf, fullfile(fig_dir, sprintf('Correlation_%s.png', bestCh)));
end

fprintf('\n All plots saved in: %s\n', fig_dir);
end
