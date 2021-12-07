function [tt1,tt2] = createTimetable(file_path)
%CREATETIMETABLE Summary of this function goes here
%   Open a .abf LFP file and export the two channels as timetables, which
%   can be further processed (pspectrum).
[d, si, ~] = abfload(file_path); 
% d : data (2 channels: d(:, 1), d(:, 2)
% si : sampling interval in µs
% ~ (not used variable, replaces h) : details of file.
s = size(d);
recording_length = s(1);
dt = si/1e6; % convert µs to s
t = (0:dt:(recording_length-1)*dt)'; % (0, dt, 2dt, ..., (recording_length-1)dt
tt1 = timetable(seconds(t), d(:, 1));
tt2 = timetable(seconds(t), d(:, 2));
end

