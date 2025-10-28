function n = count_subjects(DataMap)
% COUNT_SUBJECTS - count subjects (assumes all channels have same subject count)
chs = keys(DataMap);
if isempty(chs), n = 0; return; end
entry = DataMap(chs{1});
if isfield(entry,'group2') && ~isempty(entry.group2)
    n = numel(entry.group1) + numel(entry.group2);
else
    n = numel(entry.group1);
end
end
