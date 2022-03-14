function [belt_struct, belt_scn] = beltProcessPipeline(path_name, belt_file_name, nikon_file_name)
%BELTPROCESSPIPELINE Summary of this function goes here
% Input:
%   path_name: path to experiment folder, ends with "\".
%   belt_file_name: output txt file of labview (xy.txt; the time stamps 
%       file is xytime.txt), without '.txt' at the end
%   nikon_file_name: time stamps from NIS Elements experiment data (txt
%       file) without '.txt' at the end
% Output:
%   belt_struct: the belt data as a structure (named fields)
%   belt_scn: the belt data in scanner (Nikon) time frame.
importPackages();

if nargin == 0
    [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps();
    belt_file_name = nikon_file_name(1:end-4); %drop '_nik' ('.txt' not included), try this file name
elseif nargin == 2
    [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name);
elseif nargin == 3
    [nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name, nikon_file_name);
end

%readcaim.m
%TODO: use openImagingSession? Or at least parts of it... which do not open
%the big nikon file.
[belt, labview_time_stamps, ~, ~, ~] = openLabViewData(path_name, belt_file_name);

%TODO: here, is belt cut, or only a cut (matched to frame count of Nikon), created?
% Further analysis assumes 100 Hz... See meterPerSeocnd and possibly
% correctLength!
%FIXME: tsscn has 1 more frame at the end! Also present in original
%pipeline (readcaim.m)
[belt, tsscn] = beltMatchToNikonStamps(belt, nikon_time_stamps, labview_time_stamps); 
belt_struct = beltMatrixToStruct(belt); 
belt_struct = beltAddScannerTimeStamps(belt_struct, tsscn); %tsscn is added to belt_struct

clear belt; %only work with belt_struct from this point!
clear tsscn;

%BeltToSCN.m
belt_struct = beltCorrectArduinoArtifacts(belt_struct);
belt_struct = beltCorrectLength(belt_struct); %supply vector of zone lengths if not 500 mm per zone, 3 zones
belt_struct = beltAddRunningProperties(belt_struct);
belt_struct = beltSpeedToMeterPerSecond(belt_struct);
belt_struct = beltSmoothenPupilSize(belt_struct);
%TODO: finish smoothen pupil size function, add last step, correlate
%Ca-signals to space during periods of running as a function!
disp(length(nikon_time_stamps.data));
belt_scn = scnCreateFromBelt(belt_struct, length(nikon_time_stamps.data)); %FIXME: belt_scn is weird!

end

