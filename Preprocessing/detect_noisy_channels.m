function [badChans, chStd] = detect_noisy_channels(X, labels)
% ==============================================================
% detect_noisy_channels.m
% ==============================================================
% Detects and returns noisy EEG channels based on z-score of
% standard deviation (a simple, fast approach).
%
% INPUTS:
%   X       — EEG data [samples × channels]
%   labels  — cell array of channel labels
%
% OUTPUTS:
%   badChans — indices of noisy channels to remove
%   chStd    — standard deviation per channel
%
% ==============================================================

fprintf('Detecting noisy channels using standard deviation...\n');

% Compute channel-wise standard deviation
chStd = std(X, 0, 1);

% Convert to z-scores
zScores = (chStd - mean(chStd)) / std(chStd);

% Threshold for outlier detection
thresh = 4;   % More lenient for demo (you can lower to 3 for real data)

badChans = find(abs(zScores) > thresh);

% Print results
if isempty(badChans)
    fprintf('No noisy channels detected.\n');
else
    fprintf('Detected %d noisy channels: %s\n', ...
        numel(badChans), strjoin(labels(badChans), ', '));
end

end
