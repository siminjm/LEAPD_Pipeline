function Xf = filter_data(AllData, Filt, Fs)
% FILTER_DATA - Bandpass + 60Hz notch + L2 norm per subject vector

if iscell(AllData)
    total = numel(AllData);
else
    total = size(AllData,2);
end
Xf = cell(total,1);
[sos,g] = utils.create_filter(Fs, Filt);
wo = 60/(Fs/2); bw = wo/35; [b_notch,a_notch] = iirnotch(wo,bw);

parfor i=1:total
    X = AllData{i};
    if size(X,2)>1, X = mean(X,2); end
    X = double(X);
    X = filtfilt(sos,g,X);
    X = filtfilt(b_notch,a_notch,X);
    nrm = norm(X); if nrm>0, X = X./nrm; end
    Xf{i} = X;
end
end
