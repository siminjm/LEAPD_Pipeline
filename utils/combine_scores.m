function combo_score = combine_scores(S)
% COMBINE_SCORES - geometric mean of odds across columns
% S: NxK scores in [0,1]
odds = S ./ max(1 - S, eps);
combo_odds = geomean(odds, 2, 'omitnan');
combo_score = combo_odds ./ (1 + combo_odds);
end
