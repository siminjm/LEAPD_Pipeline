function target_vec = fetch_targets(subject_ids, labelsMap)
% FETCH_TARGETS - Match subject IDs to clinical targets.
%
% Inputs:
%   subject_ids : cell or string array of subject IDs
%   labelsMap   : containers.Map mapping ID -> target value
%
% Outputs:
%   target_vec  : numeric vector (NaN for unmatched IDs)

    target_vec = NaN(numel(subject_ids), 1);
    for i = 1:numel(subject_ids)
        sid = string(subject_ids(i));
        if isKey(labelsMap, sid)
            target_vec(i) = labelsMap(sid);
        end
    end
end
