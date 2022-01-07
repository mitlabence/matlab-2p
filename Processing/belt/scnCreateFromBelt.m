function scn = scnCreateFromBelt(belt_struct,n_nikon_frames)
%SCNCREATEFROMBELT This function creates a scanner time frame version of
%the belt_struct.
% INPUT:
%   belt_struct: belt data after processing steps
%   n_nikon_frames: if provided, the belt data is cut to match the length.
%   It is the number of frames of the .nd2 file.
% OUTPUT:
%   scn: 
%
% This is the last part of Martin's BeltToSCN.m function, with two
% differences:
%   1. This does not support NCMF customly added field "bulk" and its
%   corresponding "traceMEC" property is not used for cutting.
%   2. belt data tsscn property is left unmodified. To perform the cut on
%   this field as well, use belt_struct.tsscn = scn.tsscn.
%
% Taken from Martin's code BeltToSCN.m



%     [P,Q] = rat(length(belt_struct.tsscn)/length(belt_struct.time));
%     distancescn = resample(belt_struct.distance,P,Q);
%     runningscn = round(resample(running,P,Q));
%     roundscn = round(resample(belt_struct.round,P,Q));

speedscn = zeros(length(belt_struct.tsscn),1);
distancescn = zeros(length(belt_struct.tsscn),1);
distancePRscn = zeros(length(belt_struct.tsscn),1);
runningscn = zeros(length(belt_struct.tsscn),1);
roundscn = zeros(length(belt_struct.tsscn),1);

for i = 1:length(distancescn)
    j = abs(belt_struct.time-belt_struct.tsscn(i)) == min(abs(belt_struct.time-belt_struct.tsscn(i)));
    distancePRscn(i) = belt_struct.distancePR(find(j,1));
    distancescn(i) = belt_struct.distance(find(j,1));
    runningscn(i) = belt_struct.running(find(j,1));
    roundscn(i) = belt_struct.round(find(j,1));
    speedscn(i) = belt_struct.speed(find(j,1)); 
end
    
% %     Speed calculation using the dervative of the distance
%     for i = 2:length(distancescn)
%         speedscn(i) = (distancescn(i)-distancescn(i-1))/(belt_struct.tsscn(i)-belt_struct.tsscn(i-1));        
%     end
% %
if nargin == 2 && n_nikon_frames < length(belt_struct.tsscn) %if provided, cut to same length as Nikon frames
    int = 1:n_nikon_frames;
else %leave belt data uncut
    int = 1:length(belt_struct.tsscn);
end
%This elseif condition requires 'bulk' property which is custom and should
%not be supported.
%elseif isfield(caim,'bulk') && isfield(caim.bulk,'traceMEC') && size(caim.bulk.traceMEC,2)<length(belt_struct.tsscn) %bulk is not a CNMF field!
% %         int = length(belt_struct.tsscn)-size(caim.bulk.traceMEC,2)+1:length(belt_struct.tsscn);
%    int = 1:size(caim.bulk.traceMEC,2);      

scn.tsscn = belt_struct.tsscn(int);
scn.rounds = roundscn(int);
scn.speed = speedscn(int);
scn.distance = distancePRscn(int);
scn.totdist = distancescn(int);
scn.running = runningscn(int);

end

