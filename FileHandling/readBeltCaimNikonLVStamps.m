function [belt, caim, nikon_time_stamps, labview_time_stamps] = readBeltCaimNikonLVStamps(path_name,file_name)
%READBELTCAIMNIKONMETA First part of Martin Pofahl's readcaim function.
% Given a folder and a file name, this function tries to imply the names of
% 4 files, and open them. If they are not found, the user is asked to
% select them within a GUI. The 4 files are:
%   1. belt time stamp file: also produced by LabView; by default, with an
%       additional "time" in the file name compared to the belt data file
%       (see 2.) Used for matching with Nikon data.
%       Example: T386.021221.1105time.txt
%   2. belt data file: all the measurements made by LabView.
%       Example: T386.021221.1105.txt
%   3. nikon timestamp data: also contained in the .nd2 files, normally
%       they are extracted in NIS Viewer (or the recording software
%       for the 2-photon setup) by right-click -> Image Properties ->
%       Recorded data -> Export: Data to file. By default, this function
%       looks for the "Ca" addition to the belt data file name.
%       Example: T386.021221.1105nik.txt
%   4. CaIm NNMF results data: a .mat file. This analysis is done 
%       (currently) in the function csem.m. It uses the CaImAn functions
%       initialize_parameters, optionally manually_refine_components,
%       update_spatial_components, update_temporal_components,
%       merge_components. Normally, a suffix of "Ca" is added to the belt
%       data file name.
%       Example: T386.021221.1105Ca.mat
% Input:
%   path_name: string, path of directory where the data is located. Ends
%       with "\"!
%   file_name: file name of belt data file name without ".txt".
% Output:
%   belt:
%   caim:
%   nikon_time_stamps:
%   labview_time_stamps:

%TODO: provide suffixes to look for!
%% Load recorded belt data by LV
% Get pathanme if not given or if data is not found

if nargin == 0
    [file_name,path_name] = uigetfile('*.txt','Choose belt time stamp file');   
end
if nargin < 1 || ~exist([path_name file_name '.txt'],'file')
    [file_name,path_name] = uigetfile('*.txt','Choose belt data file',path_name);
    file_name = file_name(1:end-4);
end
disp(['Reading ' path_name file_name])

belt = importdata([path_name file_name '.txt']);

%% Load timestamps recorded by LV 

if ~exist([path_name file_name 'time.txt'],'file')
    [file_name,path_name] = uigetfile('*.txt','Choose belt time stamp file',path_name);
    labview_time_stamps = importdata([path_name file_name]);
    file_name = file_name(1:end-8);
else
    labview_time_stamps = importdata([path_name file_name 'time.txt']);
end

%% Load timestamps recorded by NIS E

if ~exist([path_name file_name 'nik.txt'],'file')
    [file_name,path_name] = uigetfile('*.txt','Choose nikon time stamp file',path_name);
    file_name = file_name(1:end-7);
end

nikon_time_stamps = importdata([path_name file_name 'nik.txt']);
nikon_time_stamps = nikon_time_stamps.data;

% Delete Artifact that sometimes occurs...
if isnan(nikon_time_stamps(end,end))  
    nikon_time_stamps = nikon_time_stamps(1:end-1,:);
end

% if realtime correction ist currupted, figure it out here
if nikon_time_stamps(2,2) == 0.1
    disp('Realtime correction is NIS Elements has been corrupted')
    nikon_time_stamps = nikon_time_stamps(:,3)*1000;
else
    nikon_time_stamps = nikon_time_stamps(:,2)*1000;
end

disp(['NIS Elements recorded frames: ' num2str(length(nikon_time_stamps))])
   

%% load Ca imaging data
disp([path_name file_name 'Ca.mat'])
if ~exist([path_name file_name 'Ca.mat'],'file')
    caim = convertComp(file_name);
else
    caim = load([path_name file_name 'Ca.mat']);
end    

if isfield(caim,'Y') && length(tsscn)<size(caim.Y,2)
    caim.Y = caim.Y(:,1:length(tsscn));
    caim.C = caim.C(:,1:length(tsscn));
    caim.S = caim.S(:,1:length(tsscn));
    caim.f = caim.f(:,1:length(tsscn));
    caim.thresh = caim.thresh(:,1:length(tsscn));
    caim.S_norm = caim.S_norm(:,1:length(tsscn));
    caim.S_bin = caim.S_bin(:,1:length(tsscn));
end

end



