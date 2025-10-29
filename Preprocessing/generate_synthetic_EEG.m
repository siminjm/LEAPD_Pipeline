function [X, Fs, labels, t] = generate_synthetic_EEG(duration_sec, nCh)
% ==============================================================
% generate_synthetic_EEG.m
% ==============================================================
% Generates realistic multi-channel EEG data with:
%   • Mixed oscillations (theta, alpha, beta)
%   • Eye blink artifacts (frontal)
%   • Slow drift (posterior)
%   • 60 Hz line noise
%   • One extremely noisy channel
%
% Suitable for testing ICA and noisy-channel detection.
%
% OUTPUTS:
%   X       — EEG signal [samples × channels]
%   Fs      — sampling rate (Hz)
%   labels  — channel labels
%   t       — time vector (seconds)
%
% Author: Simin Jamshidi (University of Iowa)
% Date: October 2025
% ==============================================================

if nargin < 1, duration_sec = 5; end
if nargin < 2, nCh = 8; end

% ----------------------------
% Basic parameters
% ----------------------------
Fs = 500;                            % Sampling rate (Hz)
t  = (0:1/Fs:duration_sec-1/Fs)';    % Time vector (seconds)

% ----------------------------
% Define neural sources
% ----------------------------
alpha = sin(2*pi*10*t);              % Alpha (10 Hz)
theta = 0.7*sin(2*pi*6*t);           % Theta (6 Hz)
beta  = 0.3*sin(2*pi*20*t);          % Beta (20 Hz)

% Three independent sources
sources = [alpha theta beta];

% ----------------------------
% Random channel mixing
% ----------------------------
mixWeights = randn(3, nCh);          % Random spatial mixing
X = sources * mixWeights;            % Mix into nCh channels

% Normalize per channel
X = X ./ max(abs(X),[],1);

% Add small Gaussian noise
X = X + 0.05*randn(size(X));

% ----------------------------
% Add typical EEG artifacts
% ----------------------------

% 1️⃣ 60 Hz line noise (global)
lineNoise = 0.1*sin(2*pi*60*t);
X = X + lineNoise;

% 2️⃣ Slow drift in posterior channels
drift = 0.2*t;
X(:, round(nCh*0.6):end) = X(:, round(nCh*0.6):end) + drift;

% 3️⃣ Eye blink artifact in frontal channels
blink = 3*exp(-((t-2.5).^2)/(2*0.12^2));   % Gaussian blink
X(:, 1:2) = X(:, 1:2) + blink;

% 4️⃣ One extremely noisy channel (simulate broken electrode)
noisyCh = nCh; % last channel
X(:, noisyCh) = X(:, noisyCh) + 6*randn(length(t),1);

% Channel labels
labels = arrayfun(@(i) sprintf('Ch%d', i), 1:nCh, 'UniformOutput', false);

% Optional quick plot
if nargout == 0
    figure('Color','w');
    plot(t, X + (0:nCh-1)*8, 'k');
    xlabel('Time (s)'); ylabel('Amplitude (µV)');
    title('Synthetic Multi-Channel EEG (demo)');
    set(gca,'YTick',[]);
end
end
