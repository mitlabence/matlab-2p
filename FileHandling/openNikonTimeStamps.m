function [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, nikon_file_name)
%OPENNIKONTIMESTAMPS Opens a txt file containing Nikon experiment data
%   (timestamps)
% Input:
%   path_name: string, path of directory where the file is located. Ends
%       with "\"!
%   nikon_file_name: file name of nikon time stamp file, without ".txt".

DISPLAY_PREFIX = "MATLAB openNikonTimeStamps: ";

if nargin == 1
        disp(strcat(DISPLAY_PREFIX, "No Nikon time stamp file name provided."));
        [nikon_file_name, path_name] = uigetfile('*.txt','Choose Nikon time stamp file', path_name);
        [~, nikon_file_name, ~] = fileparts(nikon_file_name); %cut into path, filename, .txt
else
    if nargin == 0
        disp(strcat(DISPLAY_PREFIX, "No arguments given."));
        [nikon_file_name, path_name] = uigetfile('*.txt','Choose Nikon time stamp file');
        [~, nikon_file_name, ~] = fileparts(nikon_file_name);
    elseif ~isfile(fullfile(path_name, strcat(nikon_file_name, '.txt')))
        
        disp(strcat(DISPLAY_PREFIX, "Specified file does not exist: ", fullfile(path_name, strcat(nikon_file_name, '.txt'))));
        [nikon_file_name, path_name] = uigetfile('*.txt','Choose Nikon time stamp file');
        [~, nikon_file_name, ~] = fileparts(nikon_file_name);
    end
    
end
full_path = fullfile(path_name, strcat(nikon_file_name, '.txt'));
disp(strcat(DISPLAY_PREFIX, "Reading Nikon time stamp file: ", full_path));
nikon_time_stamps = importdata(full_path);

end