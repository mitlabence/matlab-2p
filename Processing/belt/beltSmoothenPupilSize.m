function belt_struct = beltSmoothenPupilSize(belt_struct)
%BELTSMOOTHENPUPILSIZE This function performs pupil size smoothing.
% Input:
%
% Code taken from Martin's BeltToSCN.m function (Smoothing of Pupil size
% section).
%
%TODO: Not tested!

if sum(isnan(belt_struct.pupil)) > 0
    disp("BELTSMOOTHENPUPILSIZE: in belt_struct, pupil attribute contains NaN. Skipping pupil data smoothening.")
    return
end

scn.pupil = nan(length(belt_struct.tsscn),1);
scn.pupilLP = nan(length(belt_struct.tsscn),1);
scn.pupilraw = nan(length(belt_struct.tsscn),1);
scn.blink = nan(length(belt_struct.tsscn),1);

if isfield(belt_struct, 'pupil') && ~isempty(find(belt_struct.pupil,1))   
    thrs = .700;
    pupil = belt_struct.pupil;
    pupil = pupil/nanmedian(pupil);
    pupil(pupil<thrs) = nan;

    %% we want to interpolate short breaks and to delete to short onsets
    % interpolate short breaks
    win = 20; % maximum break length, if detection is unstable
    c = isnan(pupil);
    PupOff = find(diff(c)==1);
    PupOn = find(diff(c)==-1);
    if ~isempty(PupOn)
        if PupOff(1) > PupOn(1);PupOn(1) = [];end
        if length(PupOff)>length(PupOn);PupOff = PupOff(1:length(PupOn));end
        for i = 1:length(PupOff)
            if PupOn(i) - PupOff(i)<win

                int = (PupOff(i):PupOn(i)+1);
                pupil(int) = interp1([int(1) int(end)],pupil([int(1) int(end)]),int);
            end
        end
    end

    %% delete too short detection periods
    win = 1000; % minimal time window in miliseconds for stable detection
    c = isnan(pupil);
    PupOff = find(diff(c)==1);
    PupOn = find(diff(c)==-1);
    if ~c(1);PupOn =  [1; PupOn];end
    if ~isempty(PupOff) && ~isempty(PupOn)
        if PupOff(1) < PupOn(1);PupOff(1) = [];end
        PupOn = PupOn(1:length(PupOff));
        for i = 1:length(PupOn)

            if (belt_struct.time(PupOff(i)) - belt_struct.time(PupOn(i)))<win
                int = (PupOn(i):PupOff(i));
                pupil(int) = NaN;
            end
        end
    end
    %
    fs = length(belt_struct.time)/(belt_struct.time(end)/1000); 
    n = 3;       
    fc = [10]; % cutoff frequency
    Wn= fc/(fs/2);
    [b,a] = butter(n,Wn,'low');
    
    fc = [2]; % cutoff frequency
    Wn= fc/(fs/2);
    [bb,aa] = butter(n,Wn,'low');
    
    c = isnan(pupil);
    pupilLP = nan(size(pupil));
    PupOff = find(diff(c)==1);
    PupOn = find(diff(c)==-1);
    if ~c(1);PupOn =  [1; PupOn];end
    if ~isempty(PupOff)
        if PupOff(1) < PupOn(1);PupOn = [1; PupOn];end
        if PupOn(end) > PupOff(end); PupOff = [PupOff; length(pupil)];end
        for i = 1:length(PupOn)
            %%
            int = (PupOn(i)+1:PupOff(i));
            pupil(int) = filtfilt(b,a,pupil(int));
            pupilLP(int) = filtfilt(bb,aa,pupil(int));
        end
    else
        pupil = filtfilt(b,a,pupil);
    end
    belt_struct.pupilraw = belt_struct.pupil;
    belt_struct.pupilLP = pupilLP;
    belt_struct.pupil = pupil;
    
    %% Pool sizes
    if ~isempty(find(belt_struct.running,1))
        scn.pupilsize = [mean(c) mean(c(belt_struct.running==1)) mean(c(belt_struct.running==0));
                 std(c) std((c(belt_struct.running==1))) std(c(belt_struct.running==0));
                 max(c) max(c(belt_struct.running==1)) max(c(belt_struct.running==0));
                 min(c) min(c(belt_struct.running==1)) min(c(belt_struct.running==0))];
    else
        scn.pupilsize = [mean(c) NaN mean(c(belt_struct.running==0));
                 std(c) NaN std(c(belt_struct.running==0));
                 max(c) NaN max(c(belt_struct.running==0));
                 min(c) NaN min(c(belt_struct.running==0))];
    end

    %% Transfer to scanner time frame
    scn.pupil = zeros(length(belt_struct.tsscn),1);
    scn.pupilraw = zeros(length(belt_struct.tsscn),1);
    for ii = 1:length(scn.pupil)
        j = abs(belt_struct.time-belt_struct.tsscn(ii)) == min(abs(belt_struct.time-belt_struct.tsscn(ii)));
        scn.pupilraw(ii)  = belt_struct.pupilraw(find(j,1)); 
        scn.pupilLP(ii)  = belt_struct.pupilLP(find(j,1)); 
        scn.pupil(ii) = belt_struct.pupil(find(j,1));        
    end
    
    %% Locate blinking
    blink = diff(scn.pupil)<-.05;
    blink(end+1) = blink(end);
    b = find(blink);
    for ii = 1:length(b)
        blink(b(ii):b(ii)+7) = true;
    end
    blink = blink(1:length(scn.pupil));
    scn.blink = blink;
else
    belt_struct.pupil = nan(length(belt_struct.time),1);
    belt_struct.pupilLP = nan(length(belt_struct.time),1);
    belt_struct.pupilraw = nan(length(belt_struct.time),1);
    
    scn.pupil = nan(length(belt_struct.tsscn),1);
    scn.pupilLP = nan(length(belt_struct.tsscn),1);
    scn.pupilraw = nan(length(belt_struct.tsscn),1);
    scn.blink = nan(length(belt_struct.tsscn),1);
    
    scn.pupilsize = nan(4,3);
end

end

