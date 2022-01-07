function [belt_struct, belt_scn] = beltProcessPipeline(path_name)
%BELTPROCESSPIPELINE Summary of this function goes here
%   Detailed explanation goes here
importPackages();

if nargin == 0
    path_name = uigetdir("Select experiment folder");
end

[nikon_time_stamps, path_name, nikon_file_name]  = openNikonTimeStamps(path_name);

belt_file_name = nikon_file_name(1:end-4); %drop '_nik' ('.txt' not included), try this file name

%readcaim.m
%TODO: use openImagingSession? Or at least parts of it... which do not open
%the big nikon file.
[belt, labview_time_stamps, ~, ~, ~] = openLabViewData(path_name, belt_file_name);

%TODO: here, is belt cut, or only a cut (matched to frame count of Nikon), created?
% Further analysis assumes 100 Hz... See meterPerSeocnd and possibly
% correctLength!
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
belt_scn = scnCreateFromBelt(belt_struct, length(nikon_time_stamps.data)); %FIXME: belt_scn is weird!

end

