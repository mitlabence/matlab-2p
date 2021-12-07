function [t, ts1,ts2] = extractTimeSeries(file_path)
%EXTRACTTIMESERIES Extract the two channels of the .abf file for the 1-channel LFP mouse model (first channel LFP, second channel movement). Requires
%signal processing toolbox!
%   The first channel (LFP) is output raw as ts1, the second is output
%   filtered as ts2 (the movement data, which has oscillations visible)
%   Output:
%   t - series of time steps, to be used for plotting.
%   ts1 - the raw LFP data
%   ts2 - the low-pass filtered movement data
%Requirements
% Signal Processing toolbox
% abfload.m function, https://www.mathworks.com/matlabcentral/fileexchange/6190-abfload

%% Read .abf file
% d : data (2 channels: d(:, 1), d(:, 2)
% si : sampling interval in µs
% ~ (not used variable, replaces h) : details of file.
[d, si, ~] = abfload(file_path); 
s = size(d);
recording_length = s(1);
dt = si/1e6; % convert µs to s

%% Create raw output data
t = (0:dt:(recording_length-1)*dt)'; % (0, dt, 2dt, ..., (recording_length-1)dt
ts1 = d(:, 1);

%% Filter movement data
% create a filter that cuts the higher part where the spikes occur (looking
% at the spectrum). It is the same for all mice, since it is
% electronics-related.
Fpass = 0.65;
Fstop = 0.7;
Ap = 1;
Ast = 30;
filt = designfilt('lowpassfir', 'PassbandFrequency', Fpass, 'StopbandFrequency', Fstop, 'PassbandRipple', Ap, 'StopbandAttenuation', Ast);
ts2 = filtfilt(filt, d(:, 2));

end

