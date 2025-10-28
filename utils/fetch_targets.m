function y = fetch_targets(ids, labelsMap)
% FETCH_TARGETS - map string IDs to numeric targets; NaN if missing
ids = string(ids(:));
y = NaN(numel(ids),1);
for i=1:numel(ids)
    key = ids(i);
    if isKey(labelsMap, key)
        y(i) = labelsMap(key);
    else
        y(i) = NaN;
    end
end
end
