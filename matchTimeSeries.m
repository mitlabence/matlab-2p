function [outputArg1,outputArg2] = matchTimeSeries(labview_file_path, lfp_file_path)
%MATCHTIMESERIES Summary of this function goes here
%% open LFP file
%   This function should be able to match 
[~, ~, ts2] = extractTimeSeries(lfp_file_path); % use other self-written function. ts2 is the movement data, low-pass filtered.
disp(length(ts2));
%[t, ts1, ts2] is the original function

%% open LabView file
% Structure of normal LabView txt output files (not ABCtime.txt):
% columns:
% 1: cumulative distance (?)
% 2: velocity
% 3: smoother cumulative distance (?)
% 4: distance on belt (resets with new round)
% 5: reflectivity
% 6: constant
% 7: number of laps 
% 8: marks new lap (0 otherwise)
% 9: time?
% 10: time since new lap started (i.e. resets every new lap)
% the rest: constant
lv_table = readtable(labview_file_path);
velocity_data = table2array(lv_table(:, 2)); % column 2 is velocity

%% Preprocess data

ts2 = ts2 - median(ts2); % remove DC from LFP
ts2 = ts2 / max(ts2); % normalize LFP
velocity_data = velocity_data / max(velocity_data); % normalize LV data

%% plot
figure;
subplot(1, 2, 1)
plot(ts2)
title('LFP')
subplot(1, 2, 2)
plot(velocity_data)
title('LabView')
%% output
outputArg1 = ts2;
outputArg2 = velocity_data;
end

%TODO: one of the files (LFP) has way more data points. Need a way to
%downsample it!

%TODO: use xcorr() to find cross-correlation, and slide several steps to
%find maximum. ALTERNATIVE: find largest peak corresponding x, and match
%this instead!!! LFP always starts after LabView, and LFP has much more
%data points.

