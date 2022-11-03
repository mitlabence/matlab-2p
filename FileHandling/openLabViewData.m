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
%   tstamps_fname: the found LabView time stamps file name without .txt extension

%FIXME: inconsistent where ".txt" is cut off and where not (returned file
%name sometimes might contain .txt at end!)

DISPLAY_PREFIX = "MATLAB openLabViewData: ";

%% Open LabView data file
if nargin == 1
    disp(strcat(DISPLAY_PREFIX, "Only path specified: ", path_name));
    [belt_file_name,path_name] = uigetfile('*.txt','Choose belt data file', path_name);
elseif nargin == 0
    disp(strcat(DISPLAY_PREFIX, "No parameters specified."));
    [belt_file_name,path_name] = uigetfile('*.txt','Choose belt data file');
else
    fpath = fullfile(path_name, strcat(belt_file_name, '.txt'));
    disp(strcat(DISPLAY_PREFIX, "Checking if exists: ", fpath));
    if ~exist(fpath,'file')
        disp('Does not exist.');
        [belt_file_name,path_name] = uigetfile('*.txt','Choose belt data file');
    end
end
disp(strcat(DISPLAY_PREFIX, "Before fileparts: ", belt_file_name));
[~, ~, ext] = fileparts(belt_file_name);
if strcmp(ext, ".txt")  % Turns out fileparts is not so smart: cuts off anything after a period.
    [~, belt_file_name, ~] = fileparts(belt_file_name); % cut off ".txt"
end
disp(strcat(DISPLAY_PREFIX, "After fileparts: ", belt_file_name));
lv_fpath = fullfile(path_name, strcat(belt_file_name, '.txt'));
disp(strcat(DISPLAY_PREFIX, "Reading belt file: ", lv_fpath));
belt = importdata(lv_fpath);

%% Open LabView time stamp file

if ~exist(fullfile(path_name, strcat(belt_file_name, 'time.txt')),'file')
    [tstamps_fname, path_name] = uigetfile('*.txt','Choose belt time stamp file',path_name);
    disp(strcat(DISPLAY_PREFIX, "Reading belt time stamps file: ", path_name, tstamps_fname, ".txt"));
    labview_time_stamps = importdata(fullfile(path_name, tstamps_fname));
    [~, tstamps_fname, ~] = fileparts(tstamps_fname); % cut off ".txt"
else
    tstamps_fname = strcat(belt_file_name, 'time');
    disp(strcat(DISPLAY_PREFIX, "Reading belt time stamps file: ", path_name, tstamps_fname, '.txt'));
    labview_time_stamps = importdata(fullfile(path_name, strcat(tstamps_fname, '.txt')));
end


end

