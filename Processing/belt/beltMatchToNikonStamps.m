function [belt, tsscn] = beltMatchToNikonStamps(belt, nikon_time_stamps, labview_time_stamps)
%BELTMATCHTONIKONSTAMPS performs corrections for acquisition
%   software-related artifacts. The second part of Martin Pofahl's readcaim.m function,
%   with 4 changes:
%   1. odorcorrect() as first step has been removed, as presumably
%       all of our data is free of the "odor software artifact".
%   2. odorcut() as last step has been removed, for similar reasons.
%   3. Nikon imaging data is not cut (if necessary) to match LabView. In
%   most cases, this is not an issue. It is probably only needed if LabView
%   stopped recording BEFORE Nikon.
%   4. does not convert to struct. beltMatrixToStruct() does this.
%   
%   Input:
%       belt: matrix variable containing all belt read out data and the time
%             stamps.
%   Output:
%       belt: the belt matrix after corrections.
%       tsscn: belt data in the Nikon scanner time frame. Neccesary as 
%       input to beltMatrixToStruct()

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

%beltMatrixToStruct(belt, tsscn);