function [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, nikon_file_name)
%OPENNIKONTIMESTAMPS Opens a txt file containing Nikon experiment data
%   (timestamps)
% Input:
%   path_name: string, path of directory where the file is located. Ends
%       with "\"!
%   nikon_file_name: file name of nikon time stamp file, without ".txt".

if nargin == 1 %no file name provided
        [nikon_file_name, path_name] = uigetfile('*.txt','Choose Nikon time stamp file', path_name);
        nikon_file_name = nikon_file_name(1:end-4); %cut .txt
else
    if nargin == 0 || ~exist([path_name nikon_file_name '.txt'],'file')
        [nikon_file_name, path_name] = uigetfile('*.txt','Choose Nikon time stamp file');
        nikon_file_name = nikon_file_name(1:end-4); %cut .txt
    end
end
disp(['Reading Nikon time stamp file: ' path_name nikon_file_name '.txt']);
nikon_time_stamps = importdata([path_name nikon_file_name '.txt']);

end