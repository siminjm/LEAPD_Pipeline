function [best_rho, best_p, polarity, scores_aligned] = pick_polarity_and_rho(scores, y)
% PICK_POLARITY_AND_RHO - maximize |rho|; polarity=+1 keeps s, -1 flips to 1-s
mask = isfinite(scores) & isfinite(y);
if sum(mask) < 3
    best_rho = NaN; best_p = NaN; polarity = +1; scores_aligned = scores; return;
end
[sr,pr] = corr(scores(mask), y(mask), 'Type','Spearman');
sf = 1 - scores; [sf_r, sf_p] = corr(sf(mask), y(mask), 'Type','Spearman');
if abs(sf_r) > abs(sr)
    best_rho = sf_r; best_p = sf_p; polarity = -1; scores_aligned = sf;
else
    best_rho = sr; best_p = pr; polarity = +1; scores_aligned = scores;
end
end
