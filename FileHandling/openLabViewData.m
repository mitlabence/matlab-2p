function [belt, labview_time_stamps, path_name, belt_file_name, tstamps_fname] = openLabViewData(path_name,belt_file_name)
%OPENBELTDATA Opens the LabView-generated data file and time stamp file.
% Input:
%   path_name: string, path of directory where the data is located. Ends
%       with "\"!
%   belt_file_name: file name of belt data file name without ".txt".
% Output:
%   belt: the LabView belt data
%   labview_time_stamps: the LabView time stamps
%   path_name: the final path name (where the time stamp file was found)
%   belt_file_name: the found belt data file name, without .txt extension
%   tstamps_fname: the found time stamps file name without .txt extension

%FIXME: inconsistent where ".txt" is cut off and where not (returned file
%name sometimes might contain .txt at end!)

%% Open LabView data file
if nargin == 1 
    disp(['openLabViewData: Only path provided: ' path_name]);
    [belt_file_name,path_name] = uigetfile('*.txt','Choose belt data file', path_name);
    belt_file_name = belt_file_name(1:end-4); %cut off .txt
else
    disp(['openLabViewData: checking if exists: ' path_name belt_file_name '.txt']);
    if nargin == 0 || ~exist([path_name belt_file_name '.txt'],'file')
        disp('Does not exist.');
        [belt_file_name,path_name] = uigetfile('*.txt','Choose belt data file');
        belt_file_name = belt_file_name(1:end-4); %cut off .txt
    end
end

disp(['Reading belt file: ' path_name belt_file_name '.txt']);
belt = importdata([path_name belt_file_name '.txt']);

%% Open LabView time stamp file

if ~exist([path_name belt_file_name 'time.txt'],'file')
    [tstamps_fname, path_name] = uigetfile('*.txt','Choose belt time stamp file',path_name);
    disp(['Reading belt time stamps file: ' path_name tstamps_fname]);
    labview_time_stamps = importdata([path_name tstamps_fname]);
    tstamps_fname = tstamps_fname(1:end-4); %cut .txt
else
    tstamps_fname = [belt_file_name 'time'];
    disp(['Reading belt time stamps file: ' path_name tstamps_fname '.txt']);
    labview_time_stamps = importdata([path_name tstamps_fname '.txt']);
end


end

