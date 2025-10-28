function save_results(path, varstruct)
fn = fieldnames(varstruct);
for i=1:numel(fn)
    eval([fn{i} ' = varstruct.(fn{i});']); %#ok<EVLDIR>
end
save(path, fn{:});
end
