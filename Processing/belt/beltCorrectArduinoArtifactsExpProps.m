function [belt_struct, params] = beltCorrectArduinoArtifactsExpProps(belt_struct, params)
%BELTCORRECTARDUINOARTIFACTSEXPPROPS Same as beltCorrectArduinoArtifacts, 
%with additional parameter exporting functionality.
% Given a belt_struct structure of belt data,
%perform corrections for Arduino-related artifacts (outliers in speed).
%   Input:
%       belt_struct: the belt data imported as a struct (i.e. the data
%       columns named). See readcaim.m from Martin, where Belt structure is
%       created from the belt matrix.
%   Output:
%       belt_struct: the same fields as the input, with possibly modified
%       entries in speed, round, distance, distancePR fields.
%
%Code from Martin's BeltToSCN.m function (first step).
%
%TODO: make threshold values input parameters input parameters with default
%values.

% threshold for negative ardunio artifacts
thresneg = -200;
threspos = 700;
numart = find(belt_struct.speed<thresneg | belt_struct.speed>threspos);

if isprop(params, "ard_thresneg")
    disp( "beltCorrectArduinoArtifactsExpProps: ard_thresneg is overwritten!");
end
if isprop(params, "ard_threspos")
    disp( "beltCorrectArduinoArtifactsExpProps: ard_threspos is overwritten!");
end
if isprop(params, "ard_n_artifacts")
    disp( "beltCorrectArduinoArtifactsExpProps: ard_n_artifacts is overwritten!");
end
params.ard_thresneg = thresneg;
params.ard_threspos = threspos;
params.art_n_artifacts = length(numart);

for i = 1:length(numart)
    belt_struct.speed(numart(i)) = belt_struct.speed(numart(i)-1);
    belt_struct.distance(numart(i):end) = belt_struct.distance(numart(i):end)+belt_struct.distance(numart(i)-1)-belt_struct.distance(numart(i));
    currnd = belt_struct.round(numart(i));
    intend = find(belt_struct.round==currnd);intend = intend(end);
    belt_struct.distancePR(numart(i):intend) = belt_struct.distancePR(numart(i):intend)+belt_struct.distancePR(numart(i)-1)-belt_struct.distancePR(numart(i));
end

end
