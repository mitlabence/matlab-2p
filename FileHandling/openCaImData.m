function [caim, path_name, caim_file_name] = openCaImData(path_name, caim_file_name)
%OPENCAIMDATA This function opens the .mat file containing the CNMF object
%after initial processing steps. Recommended saving the object after the 
%steps in CaImAn sample: demo_script_class.m. The convention is to name this
%.mat file "xyCa.mat", where "xy.txt" is the LabView data file name. Then
%"xytime.txt" is the LabView time stamp file, "xy_nik.txt" the Nikon
%Experiment Data (time stamps) file.
% Input:
%   path_name: string, path of directory where the data is located. Ends
%       with "\"!
%   file_name: file to look for, without '.mat' extension at the end.
%TODO: NOT TESTED!


if nargin == 1 
    [caim_file_name, path_name] = uigetfile('*.txt','Choose CaImAn .mat file', path_name);
    caim_file_name = caim_file_name(1:end-4); %cut off .mat
else
    if nargin == 0 || ~exist([path_name caim_file_name '.mat'],'file')
        [caim_file_name,path_name] = uigetfile('*.mat','Choose CaImAn .mat file');
        caim_file_name = caim_file_name(1:end-4); %cut off .txt
    end
end

disp(['Reading CaImAn .mat file: ' path_name caim_file_name '.mat']);
caim = importdata([path_name caim_file_name '.mat']);

%matching Y and scanner time frame needs to be done before caim can be used
%with belt! This is the code in readcaim (appears in
%preprocessBeltCaimInplace):
%
% if isfield(caim,'Y') && length(tsscn)<size(caim.Y,2)
%    disp('Shortening following to match scanner time frame: Y, C, S, f, thresh, S_norm, S_bin');
%    disp(['reason: tsscn field (' length(tsscn) ') is shorter than Y field (' size(caim.Y,2) ').']);
%    caim.Y = caim.Y(:,1:length(tsscn));
%    caim.C = caim.C(:,1:length(tsscn));
%    caim.S = caim.S(:,1:length(tsscn));
%    caim.f = caim.f(:,1:length(tsscn));
%    caim.thresh = caim.thresh(:,1:length(tsscn));
%    caim.S_norm = caim.S_norm(:,1:length(tsscn));
%    caim.S_bin = caim.S_bin(:,1:length(tsscn));
%end
end





