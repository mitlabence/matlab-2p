function [belt,scn] = BeltToSCN(caim,belt,setlgth)
%% This function tranfers the belt data into the scanner time frame
%  On the way it corrects for some readout artifacts

% Input:
% belt          Belt data after read out
% caim          Calcium imaging data after NNMF
% setlgth       Vector with distances of stripes, if left empty length of 1500 mm and equal distances are assumed
%
% Output:
% belt          Same as input but with corrected lengths and removed artifacts
% scn           Belt data in scanner time frame
% scn.tsscn     timestamps
% scn.rounds    number of round
% scn.speed     speed of animal in m/s
% scn.distance  distance per round
% scn.totdist   total distance run in experiment 
% scn.running   logical if animal was running

%% Correct Arduino readout artifacts

% threshold for negative ardunio artifacts
thresneg = -200;
threspos = 700;
numart = find(belt.speed<thresneg | belt.speed>threspos);

for i = 1:length(numart)
    belt.speed(numart(i)) = belt.speed(numart(i)-1);
    belt.distance(numart(i):end) = belt.distance(numart(i):end)+belt.distance(numart(i)-1)-belt.distance(numart(i));
    currnd = belt.round(numart(i));
    intend = find(belt.round==currnd);intend = intend(end);
    belt.distancePR(numart(i):intend) = belt.distancePR(numart(i):intend)+belt.distancePR(numart(i)-1)-belt.distancePR(numart(i));
end

%% Normalization of belt length

stripe = find(diff(belt.stripes))+1;
rnd = find(diff(belt.round));
numstripes = max(belt.stripesPR);
% length of the belt to be normalized to
if nargin < 3 || isempty(setlgth)
    beltlgth = 1500;
    setlgth = beltlgth/numstripes:beltlgth/numstripes:beltlgth;
end
 if ~isempty(rnd)
    
    for j = belt.stripesPR(1)+1:numstripes 
        i = j-belt.stripesPR(1);
        if i ==1
            win = 1:stripe(i);
            if  belt.distancePR(win(end)) < setlgth(j)
                belt.distancePR(win) = belt.distancePR(win) + setlgth(j) - belt.distancePR(stripe(i));
            else
                belt.distancePR(win) = belt.distancePR(win)*setlgth(j)/max(belt.distancePR(win));
                belt.distance(win) = belt.distancePR(win)*setlgth(j)/max(belt.distancePR(win));
                belt.distance(win(end)+1:end) = belt.distance(win(end)+1:end) - (belt.distance(win(end)+1)-belt.distance(win(end)));
            end
        else
            win = stripe(i-1)+1:stripe(i);
            offset = belt.distancePR(win(1));                               % Offset of later steps
            belt.distancePR(win) = belt.distancePR(win) - offset;           % Thats is substracted
            fac = (setlgth(j)-setlgth(j-1))/max(belt.distancePR(win));
            belt.distancePR(win) = fac*belt.distancePR(win) + belt.distancePR(win(1)-1);
        end
        
    end

    % normalize all other rounds
    for i = 1:length(rnd)-1
        for j = 1:numstripes        
            win = stripe(i*numstripes+j-1-belt.stripesPR(1))+1:stripe(i*numstripes+j-belt.stripesPR(1));      % define the window in one round between two stripes          
            % Correction for distance per round
            if j == 1
                fac = setlgth(j)/max(belt.distancePR(win));                 % first round starts at zero, therefore only correction factor needed
                belt.distancePR(win) = fac*belt.distancePR(win);            % That is multiplied
            else
                offset = belt.distancePR(win(1));                           % Offset of later steps
                belt.distancePR(win) = belt.distancePR(win) - offset;       % Thats is substracted
                fac = (setlgth(j)-setlgth(j-1))/max(belt.distancePR(win));  % And the factor calculated for this new range
                belt.distancePR(win) = fac*belt.distancePR(win) + belt.distancePR(win(1)-1); % The last poinst of the former round is added
            end
            % Corrections for summed distance 
            offset = belt.distance(win(1));                             % Offset of later steps
            belt.distance(win) = belt.distance(win) - offset;           % Thats is substracted
            if j == 1
                fac = (setlgth(j))/max(belt.distance(win));    % And the factor calculated for this new range
            else
                fac = (setlgth(j)-setlgth(j-1))/max(belt.distance(win));    % And the factor calculated for this new range
            end
            belt.distance(win) = fac*belt.distance(win) + belt.distance(win(1)-1); % The last poinst of the former round is added
        end
    end
%%
    % normalize last round only for the summed distance
    belt.distance(rnd(end)+1:end) = belt.distance(rnd(end)+1:end)-belt.distance(rnd(end)+1)+belt.distance(rnd(end));
end

%% define periods of running

% Threshold for input signal
thres = 40;
% Width of interrunning period to be considered still running (in bins @ 100 Hz)
wind = 250;
running = zeros(length(belt.speed),1);
running(smooth(belt.speed)>thres) = 1;
stepup = find(diff(running)==1);
stepdown = find(diff(running)==-1); 
for i = 1 : length(stepdown)-1
    if length(stepup) < i+1
        
    elseif stepup(i+1)-stepdown(i)<wind
        running(stepdown(i):stepup(i+1)) = 1;
    end
end
runtime = belt.time(running == 1);

belt.running = running;
belt.runtime = runtime;

%% Transfer speed readout into m/s
convfact = 100; % factor multiplied by LabView
belt.speed = belt.speed/convfact; % 
belt.speed(1) = 0;
for i = 2:length(belt.speed)
    belt.speed(i) = belt.speed(i)/(belt.time(i)-belt.time(i-1)); %value mm/ms
%     belt.speed(i) = (belt.distance(i)-belt.distance(i-1))/(belt.time(i)-belt.time(i-1));
end

% %% Smoothing of Pupil size
% % % 
  if isfield(belt, 'pupil') && ~isempty(find(belt.pupil,1))
      thrs = 1500;
      a = belt.pupil(belt.pupil>thrs);
      b = belt.time(belt.pupil>thrs);  
      c = interp1(b,a,belt.time);
      c(isnan(c))=mean(c(~isnan(c)));
      c = smooth(c,20);
      belt.pupil = c;
     % 
      if ~isempty(find(running,1))
        scn.pupilsize = [mean(c) mean(c(running==1)) mean(c(running==0));
                 std(c) std((c(running==1))) std(c(running==0));
                  max(c) max(c(running==1)) max(c(running==0));
                 min(c) min(c(running==1)) min(c(running==0))];
     else
         scn.pupilsize = [mean(c) NaN mean(c(running==0));
                  std(c) NaN std(c(running==0));
                  max(c) NaN max(c(running==0));
                 min(c) NaN min(c(running==0))];
     end
    scn.pupil = zeros(length(belt.tsscn),1);
    for ii = 1:length(scn.pupil)
         j = abs(belt.time-belt.tsscn(ii)) == min(abs(belt.time-belt.tsscn(ii)));
        scn.pupil(ii) = belt.pupil(find(j,1));        
    end
     
 else
     scn.pupilsize = nan(4,3);
 end

%% correlate Ca-signals to space during periods of running    

%     [P,Q] = rat(length(belt.tsscn)/length(belt.time));
%     distancescn = resample(belt.distance,P,Q);
%     runningscn = round(resample(running,P,Q));
%     roundscn = round(resample(belt.round,P,Q));
    
speedscn = zeros(length(belt.tsscn),1);
distancescn = zeros(length(belt.tsscn),1);
distancePRscn = zeros(length(belt.tsscn),1);
runningscn = zeros(length(belt.tsscn),1);
roundscn = zeros(length(belt.tsscn),1);

for i = 1:length(distancescn)
    j = abs(belt.time-belt.tsscn(i)) == min(abs(belt.time-belt.tsscn(i)));
    distancePRscn(i) = belt.distancePR(find(j,1));
    distancescn(i) = belt.distance(find(j,1));
    runningscn(i) = running(find(j,1));
    roundscn(i) = belt.round(find(j,1));
    speedscn(i) = belt.speed(find(j,1)); 
end
    
% %     Speed calculation using the dervative of the distance
%     for i = 2:length(distancescn)
%         speedscn(i) = (distancescn(i)-distancescn(i-1))/(belt.tsscn(i)-belt.tsscn(i-1));        
%     end
%%

if isfield(caim,'f') && length(caim.f)<length(belt.tsscn)
%         int = length(belt.tsscn)-length(caim.f)+1:length(belt.tsscn);
    int = 1:length(caim.f);                    
elseif isfield(caim,'bulk') && isfield(caim.bulk,'traceMEC') && size(caim.bulk.traceMEC,2)<length(belt.tsscn)
%         int = length(belt.tsscn)-size(caim.bulk.traceMEC,2)+1:length(belt.tsscn);
    int = 1:size(caim.bulk.traceMEC,2);      
else
    int = 1:length(belt.tsscn);
end


scn.tsscn = belt.tsscn(int);
scn.rounds = roundscn(int);
scn.speed = speedscn(int);
scn.distance = distancePRscn(int);
scn.totdist = distancescn(int);
scn.running = runningscn(int);
belt.tsscn = belt.tsscn(int);

end