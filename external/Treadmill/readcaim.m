function [Belt,caim] = readcaim(pathname,filename)
%% Load recorded belt data by LV

if nargin == 0
    [filename,pathname] = uigetfile('*.txt','Choose belt time stamp file');   
end
if nargin < 1 || ~exist([pathname filename '.txt'],'file')
    [filename,pathname] = uigetfile('*.txt','Choose belt file',pathname);
    filename = filename(1:end-4);
end

disp(['Reading ' pathname filename])
%%

belt = importdata([pathname filename '.txt']);


%% Load timestamps recorded by LV 

if ~exist([pathname filename 'time.txt'],'file')
    [filename,pathname] = uigetfile('*.txt','Choose belt time stamp file',pathname);
    stmps = importdata([pathname filename]);
    filename = filename(1:end-8);
else
    stmps = importdata([pathname filename 'time.txt']);
end




%% Load timestamps recorded by NIS E

if ~exist([pathname filename 'nik.txt'],'file')
    [filename,pathname] = uigetfile('*.txt','Choose nikon time stamp file',pathname);
    filename = filename(1:end-7);
end

nik = importdata([pathname filename 'nik.txt']);
nik = nik.data;

% Delete Artifact that sometimes occurs...
if isnan(nik(end,end))  
    nik = nik(1:end-1,:);
end

% if realtime correction ist currupted, figure it out here
if nik(2,2) == 0.1
    disp('Realtime correction is NIS Elements has been corrupted')
    nik = nik(:,3)*1000;
else
    nik = nik(:,2)*1000;
end

disp(['NIS Elements recorded frames: ' num2str(length(nik))])
   


%% Correct for read out artifact from odor stimulation (should be corrected and can hopefully be removed soon)

belt = odorcorrect(belt);
%% Read out relevant information and compare timestamps



% Check which scanner was running and read out timestamps of scanner recorded by TF
if find(stmps(2:end,3),1)
    tmscn = stmps(stmps(:,3)~=0,3);
    disp('Scanner: galvo')
    disp(['Recorded timestamps from the microscope: ' num2str(length(tmscn))])
else
    tmscn = stmps(stmps(:,2)~=0,2);
    disp('Scanner: resonant')
    disp(['Recorded timestamps from the microscope: ' num2str(length(tmscn))])   
end
if isempty(tmscn)
    tmscn = timestampcorrect(nik,belt);
end

% Compare timestamps and data
if length(nik)>length(tmscn)
    disp(['Missed frames: ' num2str(length(nik-nik(1))-length(tmscn))])
else
    disp('No frames missed')
end


%%
if find(stmps(:,1),1)
    tmblt = stmps(stmps(:,1)~=0,1);
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


%%
% plot(belt(1:length(stmps(:,1)),9),'b')
% hold on
% plot(stmps(:,1),'r')
% hold off
% plot(belt(1:length(stmps(:,1)),9)-stmps(:,1))


%% Get belt-index interval of scanner activity

start = 1+ length(tmblt(((tmblt<tmscn(1))>0 )));
stop = start-1 + length(tmblt(((tmblt>=tmscn(1)) .* (tmblt<=tmscn(end)))>0)) ;
int = start : stop; 

%% make new timestamps for the scanner using nis elemnts time stamps
if length(nik)> length(tmscn)
    disp('NisElements timestamps used')
    tsscn =  nik-nik(1);
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

%% put belt output in struct variable for convinient handling

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

%% load Ca imging data

if ~exist([pathname filename 'Ca.mat'],'file')
    caim = convertComp(filename);
else
    load([pathname filename 'Ca.mat'])
end    

if isfield(caim,'Y') && length(tsscn)<size(caim.Y,2)
    caim.Y = caim.Y(:,1:length(tsscn));
    caim.C = caim.C(:,1:length(tsscn));
    caim.S = caim.S(:,1:length(tsscn));
    caim.f = caim.f(:,1:length(tsscn));
    caim.thresh = caim.thresh(:,1:length(tsscn));
    caim.S_norm = caim.S_norm(:,1:length(tsscn));
    caim.S_bin = caim.S_bin(:,1:length(tsscn));
end

%% Cut out stimulation period for long odor experiments

[Belt,caim] = odorcut(Belt,caim);

end

function tmscn = timestampcorrect(nik,belt)
    if find(belt(:,20),1)
        %%
        
%         disp('Realtime feedback from microscope did not work. Time will start at first pupil detection')
%         puponset = find(belt(:,20));
%         pupdiff = diff(puponset);
%         for i = 1:length(puponset)
%             if pupdiff(i) == pupdiff(i+1) == 1
%                 puponset = puponset(i);
%                 break
%             end
%         end
%         tmstrt = belt(puponset,9); %timepoint of first pupil detection in belt measurement
%         tmscn = nik + tmstrt;
        
        disp('Realtime feedback from microscope did not work. Time is alighned to last pupil detection')
        puponset = find(belt(:,20));
%         pupdiff = diff(puponset);
        pupoffset = puponset(end-1);
        tmstop = belt(pupoffset,9); %timepoint of last pupil detection in belt measurement
        tmscn = flip(tmstop - nik);
        if tmscn(1)<0
            pupoffset = puponset(1);
            tmstart = belt(pupoffset,9); %timepoint of first pupil detection in belt measurement
            tmscn = tmstart + nik;
        end
        
        figure
        plot(belt(:,9),belt(:,20))
        hold on
%         plot([tmscn(end) tmscn(end)],[0 max(belt(:,20))])
        plot([tmscn(1) tmscn(1)],[0 max(belt(:,20))],[tmscn(end) tmscn(end)],[0 max(belt(:,20))])
        hold off
        
    end
end

function belt = odorcorrect(belt)

% if isempty(find(belt(40:end-2,15:19),1))
%     return
% end

time = belt(:,9);
a = max(time);
b = find(time == a);
%%
if length(time)>b+1 && time(b +100) == 0
    belt(1:b,15:19) = belt(end-b+1:end,15:19);
    belt = belt(1:b,:);
    disp('Belt file was corrected for weird odor software artifact')   
end
%%
% plot(belt(:,9)/1000/60)
% hold on
% plot(belt(:,15)*20)
% plot(belt(:,16)*20)
% plot(belt(:,17)*30)
% plot(belt(:,18)*40)
% plot(belt(:,19)*50)
% grid on
% hold off

end

function [Belt,caim] = odorcut(Belt,caim)
    
if max(Belt.tsscn)/1000/60 < 25 || ~isfield(Belt,'odor1') || isempty(find(Belt.odor1,1))
    return
end
%%
starttime = 10; %starttime in minutes
starttime = starttime*60*1000;

startbelt = find(min(abs(Belt.time-starttime))==abs(Belt.time-starttime));
startscn = find(min(abs(Belt.tsscn-starttime))==abs(Belt.tsscn-starttime));

caim.Y = caim.Y(:,startscn:end);
caim.C = caim.C(:,startscn:end);
caim.S = caim.S(:,startscn:end);
caim.f = caim.f(:,startscn:end);
caim.thresh = caim.thresh(:,startscn:end);
caim.S_norm = caim.S_norm(:,startscn:end);
caim.S_bin = caim.S_bin(:,startscn:end);

Belt.tsscn = Belt.tsscn(startscn:end);
    Belt.tsscn = Belt.tsscn-Belt.tsscn(1);
Belt.round = Belt.round(startbelt:end);
    Belt.round = Belt.round - Belt.round(1);
Belt.speed = Belt.speed(startbelt:end);
Belt.distance = Belt.distance(startbelt:end);
    Belt.distance = Belt.distance - Belt.distance(1);
Belt.distancePR = Belt.distancePR(startbelt:end);
Belt.reflect = Belt.reflect(startbelt:end);
Belt.licking = Belt.licking(startbelt:end);
Belt.stripes = Belt.stripes(startbelt:end);
    Belt.stripes = Belt.stripes - Belt.stripes(1);
Belt.stripesPR = Belt.stripesPR(startbelt:end);
Belt.time = Belt.time(startbelt:end);
    Belt.time = Belt.time - Belt.time(1);
Belt.timePR = Belt.timePR(startbelt:end);
Belt.reward = Belt.reward(startbelt:end);
Belt.airpuff = Belt.airpuff(startbelt:end);
Belt.soundl = Belt.soundl(startbelt:end);
Belt.soundr = Belt.soundr(startbelt:end);


Belt.odor1 = Belt.odor1(startbelt:end);
Belt.odor2 = Belt.odor2(startbelt:end);
Belt.odor3 = Belt.odor3(startbelt:end);
Belt.odor4 = Belt.odor4(startbelt:end);
Belt.odor5 = Belt.odor5(startbelt:end);


if isfield(Belt,'pupil')
    Belt.pupil = Belt.pupil(startbelt:end);
end

disp(['Data has been shortened from starting minute ' num2str(starttime/1000/60) ' to a length of ' num2str(Belt.tsscn(end)/1000/60) ' min'])

end