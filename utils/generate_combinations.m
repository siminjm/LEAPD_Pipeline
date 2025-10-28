function combos = generate_combinations(best_prev, Nch, k)
% GENERATE_COMBINATIONS - progressive fixed-best using seeds from <=5-ch best
if isempty(best_prev)
    combos = nchoosek(1:Nch, k); % fallback
    return;
end
seed = [];
for t = 1:numel(best_prev)
    if isfield(best_prev(t),'indices') && ~isempty(best_prev(t).indices)
        seed = unique([seed, best_prev(t).indices]);
    end
end
seed = unique(seed);
remaining = setdiff(1:Nch, seed);
need = k - numel(seed);
if need <= 0
    combos = seed(1:k);
    return;
end
if isempty(remaining)
    combos = seed(1:min(end,k));
    return;
end
rem_combs = nchoosek(remaining, need);
combos = zeros(size(rem_combs,1), k);
for i=1:size(rem_combs,1)
    combos(i,:) = [seed, rem_combs(i,:)];
end
end
