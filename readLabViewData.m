function [outputArg1,outputArg2] = readLabViewData(inputArg1,inputArg2)
%READLABVIEWDATA Summary of this function goes here
%   Detailed explanation goes here


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
outputArg1 = inputArg1;
outputArg2 = inputArg2;
end

