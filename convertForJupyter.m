%% Description
% This script requires the LabView-generated 2 files and the Nikon
% Experiment Data txt file. From these, it does the matching of LabView to
% Nikon.
%
%

%% Clear all variables in workspace
clear; 
%% Include packages
importPackages();
%% Load data

[belt, path_name, labview_time_stamps, belt_file_name, tstamps_fname] = openLabViewData();

[nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, [belt_file_name '_nik']);

belt = correctBelt(belt, nikon_time_stamps, labview_time_stamps); 

%[belt,scn] = BeltToSCN(caim,belt); %scn is the belt in scn time frame
%scn_reduced = rmfield(scn, 'pupilsize'); %all other fields should have matching dimensions, hence writable to a single csv
%writetable(struct2table(scn_reduced), 'D:\PhD\Data\T386_MatlabTest\scn_reduced.csv');