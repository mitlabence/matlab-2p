function [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, nikon_file_name)
%OPENNIKONTIMESTAMPS Opens a txt file containing Nikon experiment data
%   (timestamps)
% Input:
%   path_name: string, path of directory where the file is located. Ends
%       with "\"!
%   nikon_file_name: file name of nikon time stamp file, without ".txt".
% Output:
%   nikon_time_stamps: now return type is table.
%       FIXME: first column is NAN! I dare not change it, as it might be
%       counted on in code further down the pipeline. Should check
%       beltProcessPipelineExpProps.m.

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
nikon_time_stamps = readtable(full_path);
% when using readtable, first column might be artificial column of NaN
% values.


% Nikon stimulation time stamp data has 6 columns. Column 5 is empty except
% for two entries: when stimulation starts, and when stimulation stops.
% The corresponding rows do NOT correspond to a frame in the nikon
% recording. Hence we must find these rows, and remove them from the stamp
% data.
% In the table, there is a NaN column as the first, so column 1 in the txt
% is column 2 in the matlab table, etc.

%check if 6th column (1st should be NaN column) exists
if any(strcmp('Var6', nikon_time_stamps.Properties.VariableNames))
    disp("MATLAB openNikonTimeStamps.m: Detected stimulation as recording type");
    % find non-empty entries in 6th column (column 5 in txt file)
    is_stimframe = zeros(height(nikon_time_stamps), 1);
    for i_frame=1:height(nikon_time_stamps)
        is_stimframe(i_frame, 1) = ~isempty(nikon_time_stamps{i_frame, 'Var6'}{1});
    end
    % is_stimframe should have 2 non-zero elements.
    % find(is_stimframe) returns a list of indices of the stimulation
    % frames (where is_stimframe is non-zero). We can remove these rows
    % from table as follows:
    nikon_time_stamps(find(is_stimframe), :) = [];

    % remove last two columns (that should be completely empty now)
    nikon_time_stamps = removevars(nikon_time_stamps, "Var6");
    nikon_time_stamps = removevars(nikon_time_stamps, "Var7");


    % Need to change Var2 from {'mm:ss.xxxx'} to sss.xxxx, double instead
    % of cell of char array
    % Also need to change Var3, Var4 from  {'xxx,yyy'} to xxx.yyy, double
    % instead of cell of char array

    % split Var2 (column of cell arrays of 'mm:ss.msms') along ':'
    min_sec_tuples = split(nikon_time_stamps{:,2}, ':');
    % add up minutes converted to seconds and seconds.milliseconds
    var1_col = nikon_time_stamps{:,1};
    var2_col = 60*cellfun(@str2num, min_sec_tuples(:,1)) + cellfun(@str2num, min_sec_tuples(:,2));
    var3_col = cellfun(@str2num, strrep(nikon_time_stamps{:,3}, ',', '.'));
    var4_col = cellfun(@str2num, strrep(nikon_time_stamps{:,4}, ',', '.'));
    var5_col = nikon_time_stamps{:,5};
    nikon_time_stamps = table(var1_col, var2_col, var3_col, var4_col, var5_col);
    end


end