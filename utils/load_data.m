function [DataMap, ChannelLoc, SubjectIDs, GroupNames] = load_data(data_file, channel_list, exclude_chans)
% LOAD_DATA - Generic loader expecting variables:
%   EEG: 1xC cell, each cell is {group1, group2?}
%   Channel_location: 1xC cellstr of names
%   Filenames (optional): struct with fields .group1 .group2 (cellstr per channel)

S = load(data_file);
ChannelLoc = [];
DataMap = containers.Map('KeyType','char','ValueType','any');
SubjectIDs = struct();
GroupNames = ["group1","group2"];

if isfield(S,'EEG') && isfield(S,'Channel_location')
    ChannelLoc = S.Channel_location;
    EEG = S.EEG;
else
    error('Expected variables EEG and Channel_location in %s', data_file);
end

if isempty(channel_list)
    channel_list = ChannelLoc;
end

% filter excluded channels
channel_list = setdiff(upper(string(channel_list)), upper(string(exclude_chans)));

for i = 1:numel(ChannelLoc)
    ch = string(ChannelLoc{i});
    if ~ismember(upper(ch), upper(channel_list)), continue; end
    cell_ch = EEG{i};
    entry = struct('group1',[],'group2',[]);
    if iscell(cell_ch) && numel(cell_ch)>=1
        entry.group1 = cell_ch{1};
    end
    if iscell(cell_ch) && numel(cell_ch)>=2
        entry.group2 = cell_ch{2};
    end
    DataMap(char(ch)) = entry;

    % Subject IDs per channel (optional if available)
    if isfield(S,'Filenames') && isstruct(S.Filenames)
        if isfield(S.Filenames,'group1')
            SubjectIDs.(char(ch)).group1 = string(S.Filenames.group1);
        else
            SubjectIDs.(char(ch)).group1 = "g1_" + (1:numel(entry.group1));
        end
        if isfield(S.Filenames,'group2')
            SubjectIDs.(char(ch)).group2 = string(S.Filenames.group2);
        else
            SubjectIDs.(char(ch)).group2 = "g2_" + (1:numel(entry.group2));
        end
    else
        % fallback: synthetic IDs
        n1 = numel(entry.group1); n2 = numel(entry.group2);
        SubjectIDs.(char(ch)).group1 = "g1_" + (1:n1);
        SubjectIDs.(char(ch)).group2 = "g2_" + (1:n2);
    end
end
end
