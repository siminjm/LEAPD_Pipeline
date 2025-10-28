function metrics = evaluate_correlation(scores, target)
% EVALUATE_CORRELATION - Spearman on complete cases
scores = scores(:); target = target(:);
mask = isfinite(scores) & isfinite(target);
if sum(mask) < 3
    metrics = struct('Rho',NaN,'PValue',NaN,'N',sum(mask));
    return;
end
[rho,p] = corr(scores(mask), target(mask), 'Type','Spearman');
metrics = struct('Rho',rho,'PValue',p,'N',sum(mask));
end
