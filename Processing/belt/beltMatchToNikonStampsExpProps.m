function [belt, tsscn, params] = beltMatchToNikonStampsExpProps(belt, nikon_time_stamps, labview_time_stamps, params)
%BELTMATCHTONIKONSTAMPSEXPPROPS same as beltMatchToNikonStamps, with an
%extra output parameter, a structure of important parameters that arose
%during the execution of the function.
%   
%   Input:
%       belt: matrix variable containing all belt read out data and the time
%             stamps.
%   Output:
%       belt: the belt matrix after corrections.
%       tsscn: belt data in the Nikon scanner time frame. Neccesary as 
%       input to beltMatrixToStruct()
%       params: struct containing which stamps were used (galvo/reso),
%       dimensions and other important parameters.

%% Compare timestamps to check for consistency

% Check which scanner was running and read out timestamps of scanner recorded by TF
%if find(labview_time_stamps(2:end,3),1)
%    tmscn = labview_time_stamps(labview_time_stamps(:,3)~=0,3);
%    disp('Scanner: galvo')
%    disp(['Recorded timestamps from the microscope: ' num2str(length(tmscn))])
%else
%    tmscn = labview_time_stamps(labview_time_stamps(:,2)~=0,2);
%    disp('Scanner: resonant')
%    disp(['Recorded timestamps from the microscope: ' num2str(length(tmscn))])   
%end
%if isempty(tmscn)
%    tmscn = timestampcorrect(nikon_time_stamps,belt);
%end

tmscn_galvo = labview_time_stamps(labview_time_stamps(:,3)~=0,3);
tmscn_reso = labview_time_stamps(labview_time_stamps(:,2)~=0,2);
disp(['Size of galvo: ' num2str(size(tmscn_galvo)) '; reso: ' num2str(size(tmscn_reso))]);
if any(size(tmscn_galvo) > size(tmscn_reso))
    tmscn = tmscn_galvo;
    disp("Galvo is longer, choosing it");
else
    disp("Reso is longer, choosing it");
    tmscn = tmscn_reso;
end

if isprop(params, "lvstamps_len")
    disp( "beltMatchToNikonStampsExpProps: lvstamps_len is overwritten!");
end
params.lvstamps_len = length(tmscn);


% Compare timestamps and data
if isprop(params, "missed_frames")
    disp( "beltMatchToNikonStampsExpProps: missed_frames is overwritten!");
end

n_missed_frames = length(nikon_time_stamps.data)-length(tmscn);
params.missed_frames = n_missed_frames; 
if n_missed_frames > 0
    disp(['Missed frames: ' num2str(n_missed_frames)]);
else 
    if n_missed_frames < 0 % more labview stamps, length(tmscn) larger
        disp(['Nikon time stamps less than LabView-recorded time stamps! ' num2str(length(nikon_time_stamps.data)) ' (Nikon) vs ' num2str(length(tmscn)) ' (LabView).']);
    else
        disp('No frames missed');
    end
end


if find(labview_time_stamps(:,1),1)
    tmblt = labview_time_stamps(labview_time_stamps(:,1)~=0,1);
    n_missed_cycles = length(belt(:,9))-length(tmblt);

    if isprop(params, "missed_cycles")
        disp( "beltMatchToNikonStampsExpProps: missed_cycles is overwritten!");
    end
    params.missed_cycles = n_missed_cycles;
    if n_missed_cycles == 0
        disp('No cycles missed')
    else
        disp(['Missed belt cycles: ' num2str(n_missed_cycles)])
        if abs(n_missed_cycles)>10  % TODO: why 10 as threshold?
            disp('Timestamps filled with the belt stamps')
            tmblt = [tmblt(1); belt(:,9)+tmblt(2)-belt(1,9)];
        end
    end
else
    if isprop(params, "missed_cycles")
        disp( "beltMatchToNikonStampsExpProps: missed_cycles is overwritten!");
    end
    params.missed_cycles = NaN;
    disp('Belt time stamps were not recorded');
    tmblt = belt(:,9);    
end

%% Get belt-index interval of scanner activity

start = 1+ length(tmblt(((tmblt<tmscn(1))>0 )));
stop = start-1 + length(tmblt(((tmblt>=tmscn(1)) .* (tmblt<=tmscn(end)))>0)) ;
int = start : stop;

if isprop(params, "i_belt_start")
    disp( "beltMatchToNikonStampsExpProps: i_belt_start is overwritten!");
end
if isprop(params, "i_belt_stop")
    disp( "beltMatchToNikonStampsExpProps: i_belt_stop is overwritten!");
end
params.i_belt_start = start;
params.i_belt_stop = stop; % inclusive stop!

%% make new timestamps for the scanner using nis elements time stamps
if length(nikon_time_stamps.data)> length(tmscn)
    disp('NisElements timestamps used')
    tsscn =  nikon_time_stamps.data-nikon_time_stamps.data(1,:);
    tsscn = tsscn(:,2)*1000;

    if isprop(params, "used_tstamps")
        disp( "beltMatchToNikonStampsExpProps: used_tstamps is overwritten!");
    end
    params.used_tstamps = "Nikon";
else %TODO: ERROR for the case more LabView stamps than Nikon stamps
    disp('LabView timestamps used')
    tsscn = tmscn - tmscn(1);

    if isprop(params, "used_tstamps")
        disp( "beltMatchToNikonStampsExpProps: used_tstamps is overwritten!");
    end
    params.used_tstamps = "LabView";
end
if(isnan(tsscn(1,1))) %remove all-NaN column
    tsscn(:, 1) = [];
end

if isprop(params, "len_tsscn")
    disp( "beltMatchToNikonStampsExpProps: len_tsscn is overwritten!");
end
params.len_tsscn = length(tsscn);

% check if timestamp vector is unique
if length(tsscn) ~= length(unique(tsscn))
    if isprop(params, "timestamps_were_duplicate")
        disp( "beltMatchToNikonStampsExpProps: timestamps_were_duplicate is overwritten!");
    end
     params.timestamps_were_duplicate = true;

     a = diff(tsscn) > 0;
     tsscnnew = interp1(find(a),tsscn(a),1:length(tsscn));
     tsscnnew(end) = tsscn(end);
     tsscn = tsscnnew'; %transpose!
     disp('Timestamps are made unique');
end


if isprop(params, "timestamps_were_duplicate")
    disp( "beltMatchToNikonStampsExpProps: timestamps_were_duplicate is overwritten!");
end
if isprop(params, "len_tsscn_new")
    disp( "beltMatchToNikonStampsExpProps: len_tsscn_new is overwritten!");
end

params.timestamps_were_duplicate = false;
params.len_tsscn_new = length(tsscn);


%TODO: what do these do? (Documentation)
tsblt = tmblt(int)-tmblt(start);
belt = belt(int,:);
belt(:,9) = tsblt;

movie_length = tsscn(end)/60/1000;
freq = length(tsscn)/tsscn(end)*1000;


if isprop(params, "movie_length_min")
    disp( "beltMatchToNikonStampsExpProps: movie_length_min is overwritten!");
end
if isprop(params, "frequency_estimated")
    disp( "beltMatchToNikonStampsExpProps: frequency_estimated is overwritten!");
end

params.movie_length_min = movie_length;
params.frequency_estimated = freq;

disp(['Movie length: ' num2str(movie_length) ' min'])
disp(['Frequency: ' num2str(freq) ' Hz'])

%beltMatrixToStruct(belt, tsscn);