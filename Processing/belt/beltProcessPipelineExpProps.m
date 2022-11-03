function [belt_struct, belt_scn, processing_params] = beltProcessPipelineExpProps(path_name, belt_file_name, nikon_file_name)
%BELTPROCESSPIPELINEEXPPROPS Same as beltProcessPipeline, with an
%additional output parameter containing important parameters.
% Input:
%   path_name: path to experiment folder, optionally ends with separator
%   ("\" or "/")
%   belt_file_name: output txt file of labview (xy.txt; the time stamps 
%       file is xytime.txt), without '.txt' at the end
%   nikon_file_name: time stamps from NIS Elements experiment data (txt
%       file) without '.txt' at the end
% Output:
%   belt_struct: the belt data as a structure (named fields)
%   belt_scn: the belt data in scanner (Nikon) time frame.
%   processing_params: struct containing missed frames, missed cycles, and
%   other important parameters that arose during the processing pipeline.
%   %TODO: document processing_params


%FIXME: if nikon time stamps file is not named xy_nik.txt but xynik.txt
%(which happens often), this function will not work properly!
importPackages();

processing_params = struct;

DISPLAY_PREFIX = "MATLAB beltProcessPipelineExpProps: "; % put this in beginning of all disp() to mark which function displays message

if nargin == 0
    % Ask user to select files
    [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps();
    belt_file_name = nikon_file_name(1:end-4); %drop '_nik' ('.txt' not included), try this file name
elseif nargin == 1
    % Only path is given, ask user to select files
    [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name);
    %try to infer belt file name from nikon file name
    if endsWith(nikon_file_name, "_nik")
        belt_file_name = nikon_file_name(1:end-4); %drop '_nik' ('.txt' not included), try this file name
    elseif endsWith(nikon_file_name, "nik")
        belt_file_name = nikon_file_name(1:end-3); %drop 'nik', try this file name
    else
        % belt file name cannot be derived from nikon time stamp minus _nik
        % or nik at the end
        
        % TODO: should handle this error with a prompt to specify belt file instead of error!
        error(strcat(DISPLAY_PREFIX, "1 argument provided, given nikon stamp file does not end with _nik or nik: ", nikon_file_name))
    end
elseif nargin == 2
    %path and belt file given, need to figure out nikon file name
    %try to infer nikon file name:
    if isfile(fullfile(path_name, strcat(belt_file_name, "_nik.txt")))
        nikon_file_name = strcat(belt_file_name, "_nik");
        disp(strcat(DISPLAY_PREFIX, "Trying to open inferred nikon stamp file name: ", nikon_file_name, ".txt"));
        [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, nikon_file_name);
    elseif isfile(fullfile(path_name, strcat(belt_file_name, "nik.txt")))
        nikon_file_name = strcat(belt_file_name, "nik");
        disp(strcat(DISPLAY_PREFIX, "Trying to open inferred nikon stamp file name: ", nikon_file_name, ".txt"));
        [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, nikon_file_name);
    else
        disp("path and belt file name provided; failed to infer Nikon time stamp file.");
        [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name);
    end
elseif nargin == 3
    disp(strcat(DISPLAY_PREFIX, "Nikon file name specified; trying to open: ", nikon_file_name, ".txt"));
    [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, nikon_file_name);
end

processing_params.path_name = path_name;
processing_params.belt_file_name = belt_file_name;
processing_params.nikon_file_name = nikon_file_name;

%readcaim.m
%TODO: use openImagingSession? Or at least parts of it... which do not open
%the big nikon file.
[belt, labview_time_stamps, ~, ~, tstamps_file_name] = openLabViewData(path_name, belt_file_name);

processing_params.tstamps_file_name = tstamps_file_name;

%TODO: here, is belt cut, or only a cut (matched to frame count of Nikon), created?
% Further analysis assumes 100 Hz... See meterPerSeocnd and possibly
% correctLength!
%FIXME: tsscn has 1 more frame at the end! Also present in original
%pipeline (readcaim.m)
[belt, tsscn, processing_params] = beltMatchToNikonStampsExpProps(belt, nikon_time_stamps, labview_time_stamps, processing_params);

belt_struct = beltMatrixToStruct(belt); 
belt_struct = beltAddScannerTimeStamps(belt_struct, tsscn); %tsscn is added to belt_struct

clear belt; %only work with belt_struct from this point!
clear tsscn;

%BeltToSCN.m
[belt_struct, processing_params] = beltCorrectArduinoArtifactsExpProps(belt_struct, processing_params);
[belt_struct, processing_params] = beltCorrectLengthExpProps(belt_struct, processing_params); %supply vector of zone lengths if not 500 mm per zone, 3 zones
[belt_struct, processing_params]  = beltAddRunningPropertiesExpProps(belt_struct, processing_params);
belt_struct = beltSpeedToMeterPerSecond(belt_struct);
belt_struct = beltSmoothenPupilSize(belt_struct);
%TODO: finish smoothen pupil size function, add last step, correlate
%Ca-signals to space during periods of running as a function!
disp(length(nikon_time_stamps.data));
belt_scn = scnCreateFromBelt(belt_struct, length(nikon_time_stamps.data)); %FIXME: belt_scn is weird!

end

