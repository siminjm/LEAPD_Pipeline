function T = read_labels_table(filepath)
% READ_LABELS_TABLE - Standardize the clinical label file
%
% Expects an Excel or CSV file with at least:
%   ID, and either Target or one of {Days, MoCA, UPDRS}.
%
% Returns:
%   T table with columns ID and Target.

    if isempty(filepath) || ~isfile(filepath)
        warning('No labels file found or path invalid: %s', filepath);
        T = table();
        return;
    end

    T = readtable(filepath, 'VariableNamingRule','preserve');

    % Standardize column names
    if ~ismember('ID', T.Properties.VariableNames)
        error('Clinical label file must contain an "ID" column.');
    end

    if ~ismember('Target', T.Properties.VariableNames)
        possible = intersect(["Days","MoCA","UPDRS","Score","Value"], ...
                             string(T.Properties.VariableNames));
        if isempty(possible)
            error('Clinical label file must contain a "Target" column or recognizable alternative.');
        end
        T.Target = T.(possible(1)); % rename automatically
    end
end
