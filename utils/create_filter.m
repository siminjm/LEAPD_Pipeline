function [sos,g] = create_filter(Fs, Filt)
% CREATE_FILTER - Butterworth bandpass defined by Filt rows [0 f1 0; f2 inf 0]
start_f = Filt(1,2);
if (size(Filt,1)==1) || (Filt(2,3)==1)
    stop_f = Fs/2 - 1e-3;
else
    stop_f = Filt(2,1);
end
Order = 6; Fn = Fs/2; Wp = [start_f stop_f]/Fn;
[z,p,k] = butter(Order, Wp, 'bandpass');
[sos,g] = zp2sos(z,p,k);
end
