function [Belt, caim] = correctBeltCaim(belt, caim, nikon_time_stamps, labview_time_stamps)
%preprocessBeltCaimInplace perform corrections for acquisition 
%   software-related artifacts. The second part of Martin Pofahl's readcaim.m function,
%   with two changes:
%   1. odorcorrect() as first step has been removed, as presumably
%       all of our data is free of the "odor software artifact".
%   2. odorcut() as last step has been removed, for similar reasons.
%   
%   Input:
%       belt: struct variable containing all belt read out data and the time
%             stamps
%       caim: output caim of csem. Unfortunately, at this point, the object
%               type cannot be determined.

%% Compare timestamps to check for consistency

% Check which scanner was running and read out timestamps of scanner recorded by TF
if find(labview_time_stamps(2:end,3),1)
    tmscn = labview_time_stamps(labview_time_stamps(:,3)~=0,3);
    disp('Scanner: galvo')
    disp(['Recorded timestamps from the microscope: ' num2str(length(tmscn))])
else
    tmscn = labview_time_stamps(labview_time_stamps(:,2)~=0,2);
    disp('Scanner: resonant')
    disp(['Recorded timestamps from the microscope: ' num2str(length(tmscn))])   
end
if isempty(tmscn)
    tmscn = timestampcorrect(nikon_time_stamps,belt);
end

% Compare timestamps and data
if length(nikon_time_stamps)>length(tmscn)
    disp(['Missed frames: ' num2str(length(nikon_time_stamps-nikon_time_stamps(1))-length(tmscn))])
else
    disp('No frames missed')
end

if find(labview_time_stamps(:,1),1)
    tmblt = labview_time_stamps(labview_time_stamps(:,1)~=0,1);
    if length(tmblt) == length(belt(:,9))       
        disp('No cycles missed')
    else
        disp(['Missed belt cycles: ' num2str(length(belt(:,9))-length(tmblt))])
        if abs(length(belt(:,9))-length(tmblt))>10
            disp('Timestamps filled with the belt stamps')
            tmblt = [tmblt(1); belt(:,9)+tmblt(2)-belt(1,9)];
        end
    end
else
    disp('Belt time stamps were not recorded')
    tmblt = belt(:,9);    
end

%% Get belt-index interval of scanner activity

start = 1+ length(tmblt(((tmblt<tmscn(1))>0 )));
stop = start-1 + length(tmblt(((tmblt>=tmscn(1)) .* (tmblt<=tmscn(end)))>0)) ;
int = start : stop; 

%% make new timestamps for the scanner using nis elements time stamps
if length(nikon_time_stamps)> length(tmscn)
    disp('NisElements timestamps used')
    tsscn =  nikon_time_stamps-nikon_time_stamps(1);
else
    disp('LabView timestamps used')
    tsscn = tmscn - tmscn(1);
end
% check if timestamp vector is unique
if length(tsscn) ~= length(unique(tsscn))
     a = diff(tsscn) > 0;
     tsscnnew = interp1(find(a),tsscn(a),1:length(tsscn));
     tsscnnew(end) = tsscn(end);
     tsscn = tsscnnew';
     disp('Timestamps are made unique')
end
tsblt = tmblt(int)-tmblt(start);
belt = belt(int,:);
belt(:,9) = tsblt;

disp(['Movie length: ' num2str(tsscn(end)/60/1000) ' min'])
disp(['Frequency: ' num2str(length(tsscn)/tsscn(end)*1000) ' Hz'])

%% create belt output in struct variable for convinient handling

Belt = struct;

Belt.tsscn = tsscn;
Belt.round = belt(:,1) - belt(1,1);
Belt.speed = belt(:,2);
Belt.distance = belt(:,3) - belt(1,3);
Belt.distancePR = belt(:,4);
Belt.reflect = belt(:,5);
Belt.licking = belt(:,6);
Belt.stripes = belt(:,7) - belt(1,7);
Belt.stripesPR = belt(:,8);
Belt.time = belt(:,9);
Belt.timePR = belt(:,10);
Belt.reward = belt(:,11);
Belt.airpuff = belt(:,12);
Belt.soundl = belt(:,13);
Belt.soundr = belt(:,14);

if size(belt,2)>14
    Belt.odor1 = belt(:,15);
    Belt.odor2 = belt(:,16);
    Belt.odor3 = belt(:,17);
    Belt.odor4 = belt(:,18);
    Belt.odor5 = belt(:,19);
end

if size(belt,2)>19
    Belt.pupil = belt(:,20);
end
%% Cut nikon imaging data if necessary
%TODO: put this in separate function!
nikonCutToLabView(Belt, caim);

end