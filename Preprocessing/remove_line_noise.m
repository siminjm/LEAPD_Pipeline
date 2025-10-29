function X_clean = remove_line_noise(X, Fs, baseHz, Q)
%REMOVE_LINE_NOISE  Notch-filter 60 Hz (and harmonics) with zero-phase.
%   X may be [samples x channels] or [channels x samples]
%   Fs = sampling rate (Hz)
%   baseHz = base line frequency (default 60)
%   Q = quality factor (default 35; higher = narrower notch)
%
%   Example:
%   X_clean = remove_line_noise(X, 500);              % 60 Hz, 120, 180...
%   X_clean = remove_line_noise(X, 512, 60, 45);      % slightly wider notch

    if nargin < 3 || isempty(baseHz), baseHz = 60; end
    if nargin < 4 || isempty(Q), Q = 35; end

    % Orient to [samples x channels]
    flipped = false;
    if size(X,1) < size(X,2)
        X = X.'; flipped = true;
    end
    X_clean = X;

    nyq = Fs/2;
    freqs = baseHz:baseHz:nyq-1e-9;   % all harmonics strictly below Nyquist

    for f0 = freqs
        w0 = f0/(Fs/2);           % normalize to [0,1] (1 -> Nyquist)
        bw = w0 / Q;              % iirnotch uses bw = w0/Q
        [b,a] = iirnotch(w0, bw);
        % zero-phase (no group delay)
        X_clean = filtfilt(b, a, X_clean);
    end

    if flipped
        X_clean = X_clean.';
    end
end
