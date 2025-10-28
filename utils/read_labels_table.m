function T = read_labels_table(labels_file)
% READ_LABELS_TABLE - load ID/Target from .xlsx or .csv (robust to header names)
[~,~,ext] = fileparts(labels_file);
switch lower(ext)
    case {'.xlsx','.xls'}
        T = readtable(labels_file, 'FileType','spreadsheet');
    case '.csv'
        T = readtable(labels_file, 'FileType','text');
    otherwise
        error('Unsupported label file: %s', labels_file);
end

% Try to standardize column names
vars = lower(string(T.Properties.VariableNames));
idIdx = find(ismember(vars, ["id","subject","subjectid","sid","name"]),1);
tgtIdx = find(ismember(vars, ["target","label","score","y","value"]),1);
if isempty(idIdx) || isempty(tgtIdx)
    error('Labels table must contain ID and Target columns (or recognizable aliases).');
end
T = T(:, [idIdx tgtIdx]);
T.Properties.VariableNames = {'ID','Target'};
end
