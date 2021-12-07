function scn = stimcor(belt,caim,scn,showres)
%%
if isfield(caim,'C')
    S = caim.S_bin;
    C = caim.C;
else
    S = [];
    C = [];
end
if nargin < 4
    showres = 1;
end

%% reward analysis
if length(find(diff(belt.reward)))>2
    wind = [-2 8];                                     % size of window to look for response
    thresh = 5;                                        % minimal number of detected events around stimulus
    respwind = [0 1];                                  % size of window to look for response
    scn.reward = findstim(belt.reward,caim,wind,thresh,respwind,S,C);
    if showres;disp(['Animal was rewarded ' num2str(length(scn.reward.timepnt)/2) ' times']);end
end

%% tone left analysis
if length(find(diff(belt.soundl)))>2
    wind = [-2 8];                                     % size of window to look for response
    thresh = 5;                                        % minimal number of detected events around stimulus
    respwind = [0 1];                                  % size of window to look for response
    scn.soundl = findstim(belt.soundl,caim,wind,thresh,respwind,S,C);  
%     stimplot(scn.soundl)    
%     set(gcf,'name','Tone Left')
    if showres;disp(['Auditory cues on left side were presented ' num2str(length(scn.soundl.timepnt)/2) ' times']);end
end

%% tone right analysis
if length(find(diff(belt.soundr)))>2
    wind = [-2 8];                                     % size of window to look for response
    thresh = 5;                                        % minimal number of detected events around stimulus
    respwind = [0 1];                                  % size of window to look for response
    scn.soundr = findstim(belt.soundr,caim,wind,thresh,respwind,S,C); 
%     stimplot(scn.soundr)    
%     set(gcf,'name','Tone Right')
    if showres;disp(['Auditory cues on right side were presented ' num2str(length(scn.soundr.timepnt)/2) ' times']);end
end

%% network response analysis
if isfield(caim,'network') && ~isempty(caim.network)
    wind = [-2 8];                                  % size of window to look for response
    thresh = 1;                                      % minimal number of detected events around stimulus
    respwind = [-.1 .3];                                 % size of window to look for response
    scn.network = findstim(caim.network.netev,caim,wind,thresh,respwind,S,C); 
%     stimplot(scn.network)
%     set(gcf,'name','network')  
end

%% airpuff analysis
if length(find(diff(belt.airpuff)))>2
    wind = [-2 8];                                  % size of window to look for response
    thresh = 4;                                      % minimal number of detected events around stimulus
    respwind = [0 1];                                 % size of window to look for response
    scn.airpuff = findstim(belt.airpuff,caim,wind,thresh,respwind,S,C);
    scn.airpuff = stimvect(scn.airpuff);
    % correct stimuli cell ID
    if isfield(scn,'pcell')
        stim = scn.airpuff;
        psig = pcell(stim.cellID(:,1),4) <= .05;
        pexc = pcell(stim.cellID(:,1),4) > .05;

        stim.nr.cellID = [stim.nr.cellID;stim.cellID(pexc,:)];
        stim.nr.resp = cat(1,stim.nr.resp, stim.resp(pexc,:,:));
        stim.nr.respCa = cat(1,stim.nr.respCa, stim.respCa(pexc,:,:));
        stim.nr.maxresppost = cat(1,stim.nr.maxresppost, stim.maxresppost(pexc,:));
        stim.nr.maxresppre = cat(1,stim.nr.maxresppre, stim.maxresppre(pexc,:));

        stim.cellID = stim.cellID(psig,:);
        stim.resp = stim.resp(psig,:,:);
        stim.respCa = stim.respCa(psig,:,:);
        stim.maxresppost = stim.maxresppost(psig,:);
        stim.maxresppre = stim.maxresppre(psig,:);
        stim.stimpol = stim.stimpol(psig,:);
    end
%     stimplot(scn.airpuff)    
%     set(gcf,'name','Airpuff')
    if showres;disp(['Airpuffs were presented ' num2str(length(scn.airpuff.timepnt)/2) ' times']);end
end

%% airpuff with running onset analysis
% if ~isempty(find(belt.airpuff,1))
%     wind = -2:8;                                     % size of window to look for response
%     thresh = 5;                                      % minimal number of detected events around stimulus
%     respwind = 0:1 ;                                 % size of window to look for response
%     
%     a = diff(belt.running);
%     a(end+1) = a(end);
%     a(a==-1) = 0;
%     win = 1:200;
%     a(win) = 0;
%     a(find(a)-win) = 1;
%     b = zeros(size(a));
%     b(belt.airpuff == 1 & a == 1) = 1;
%     
%     scn.airpuffrun = findstim(b,caim,wind,thresh,respwind,S,C);
%     b = zeros(size(a));
%     b(belt.airpuff == 1 & a == 0) = 1;
%     scn.airpuffnorun = findstim(b,caim,wind,thresh,respwind,S,C);
%     disp(['Airpuffs triggered running ' num2str(length(scn.airpuffrun.timepnt)/2) ' times'])   
% end

%% running onset analysis
if length(find(diff(belt.running)))>2
    wind = [-2 8];                                     % size of window to plot around response
    thresh = 3;                                      % minimal number of detected events around stimulus
    respwind = [-.5 .5];                               % size of window to look for response
    scn.runonset = findstim(belt.running,caim,wind,thresh,respwind,S,C); 
    if showres;disp(['Mouse started running ' num2str(length(scn.runonset.timepnt)/2) ' times']);end
%     stimplot(scn.runonset) 
%     set(gcf,'name','Run onset')
end

%% place coding analysis
if isfield(scn,'spcon') && ~isempty(find(scn.spcon(4,:),1))
    wind = [-2 8];                                     % size of window to plot around response
    thresh = ceil(.2*max(scn.rounds));                 % minimal number of detected events around stimulus
    respwind = [0 .2];                                 % size of window to look for response
    scn.spcresp = findstim(scn.spcon(4,:),caim,wind,thresh,respwind,S,C); 
    if showres;disp(['Place related firing occured ' num2str(length(scn.spcresp.timepnt)/2) ' times']);end
%     stimplot(scn.spcresp) 
%     set(gcf,'name','Place code')
end

%% belt stripe analysis
if length(find(diff(belt.stripesPR)))>2 && max(belt.stripesPR)>1
    wind = [-2 8];                                    % size of window to look for response
    thresh = 4;                                     % minimal number of detected events around stimulus
    respwind = [0 1]; 
    stripes = belt.stripesPR;
    stripes(2:end) = diff(stripes);
    stripes(belt.stripesPR ~=2) = 0;
    scn.stripes = findstim(stripes,caim,wind,thresh,respwind,S,C); 
    if showres;disp(['Passed stripe 2 ' num2str(length(scn.stripes.timepnt)/2) ' times']);end
end


function Resp = findstim(stim,caim,wind,thresh,respwind,S,C)
%   inputs: 
%       stim        stimuli in belt framework
%       caim        Ca imaging Data
%       S           binary Ca onsets
%       C           deconvolved Ca-traces
%       wind        window for plotting porpuses
%       thresh      number of minimal responses
%       respwind    window where to look for responses around stimuli

    Fs = 1000/scn.tsscn(end)*length(scn.tsscn);
    wind = round(wind(1)*Fs):round(wind(end)*Fs);
    respwind = round(respwind(1)*Fs):round(respwind(end)*Fs);
    
    %% find stimulation timepoints in scanner timestamps 
    
    if length(stim)>length(scn.tsscn)
        stretch = ceil(length(belt.time)/length(scn.tsscn));
        stim([1:(abs(wind(end))*stretch+stretch) end-(abs(wind(end))*stretch):end]) = 0;  
        timesbelt = find(diff(stim));                       % Find stimulation points in belt timeframe   
        timesscn = zeros(length(timesbelt),1);              % Translated timepoints in scanner timeframe
        stimscn = zeros(length(scn.tsscn),1);              % Vector of stim activity
        exclude = ones(1,length(timesscn));

        for i = 1:2:length(timesbelt)
            [~,j]  = min(abs(scn.tsscn-belt.time(timesbelt(i))));
            [~,jj] = min(abs(scn.tsscn-belt.time(timesbelt(i+1))));
            if j > abs(wind(1)) && jj < length(scn.tsscn)-abs(wind(end))              
                stimscn(j:jj) = 1;
            else 
                exclude(i) = 0;exclude(i+1) = 0;
            end
        end
        timesscn = timesscn(exclude==1);
        numstim = sum(diff(stimscn)==1);          % Number of stimulations
    else
        stim([1:wind(end) end-wind(end):end]) = 0;
        timesscn = scn.tsscn(diff(stim)==1|diff(stim)==-1);
        numstim = length(timesscn)/2;
        stimscn = stim;
    end
    timesscn = scn.tsscn(diff(stimscn)==1|diff(stimscn)==-1);
    pretime = timesscn(1:2:end);
    pretime = [0; diff(pretime)/1000];
    
    Resp.stimon = stimscn;
    Resp.timepnt = timesscn;
    Resp.pretime = pretime;
    %% Read out positions of stimulation in binned belt
    
    spcbin = 20; % size of bin in mm
%     numbin = round(max(belt.distancePR)/spcbin);
    numbin = round(1500/spcbin);
    stimpos = scn.distance(diff(stimscn)==1);
    stimbin = round(stimpos./spcbin);
    
    if max(belt.round)>1           
%         
%         stimplc = zeros(1,numbin);
%         for i = 1:numbin
% %             stimplc(i) = sum(stimpos > (i-1)*spcbin+1 & stimpos < i*spcbin);              
%             wind(i) = sum(stimpos==i);
%         end
%         stimplc = histcounts(round(stimpos./spcbin),1:101);
        
        stimplc = histcounts(stimbin,1:numbin+1);
        if length(stimplc)>100;stimplc = stimplc(1:100);end       
    else
        stimplc = NaN(1,numbin);
        
    end
    Resp.stimplc = stimplc;
    
    %% Extract responses in time window around stim points
    
    if isfield(caim,'C')    %  extract Ca-traces for every point and every cell
        Data = C;
        DecData = S;
        numcells = length(Data(:,1));                       % Number of cells
        resp = zeros(numcells,numstim,length(wind));        % binary Data to look for responses
        respCa = zeros(numcells,numstim,length(wind));      % traces of Data to look for responses        
        respPlCl = zeros(numcells);                         % look only at response of place coding cells
        respSpec = ones(numcells,numstim);                  % look on responders, when they are not responding
        cellplc = zeros(numcells,numstim);                  % positions of responded stimuli
        numresp = zeros(numcells,2);                        % number of responses
        maxAmppost = zeros(numcells,numstim);               % max response after stimulus
        maxAmppre = zeros(numcells,numstim);                % max response before stimulus               
        cellID = zeros(numcells,1);
        cellIDnr = zeros(numcells,1);
        respwind = find(wind == respwind(1)):find(wind == respwind(1))+size(respwind,2)-1;
        
        for i = 1:numcells 
            for j = 1:numstim               % Go through all onsets of stimuli to get the responses
                resp(i,j,:)   = DecData(i,find(scn.tsscn == timesscn(2*j-1))+wind);
                respCa(i,j,:) = Data(i,find(scn.tsscn == timesscn(2*j-1))+wind);              
                if find(resp(i,j,respwind),1)
                    numresp(i,1)  = numresp(i,1) + sum(resp(i,j,respwind));
                    numresp(i,2)  = numresp(i,2) + 1;
                    cellplc(i,j) = stimbin(j);
                    respSpec(i,j) = 0;
                end
                maxAmppost(i,j) = max(respCa(i,j,find(wind==0,1):end));
                maxAmppre(i,j)  = max(respCa(i,j,1:find(wind==0,1)));
                
            end           
            
            if numresp(i,1) > abs(thresh)
                cellID(i) = i;
            else
                cellIDnr(i) = i; 
                respSpec(i,:) = 0;
            end
            
            if isfield(scn,'cellID') && ~isempty(find(scn.cellID(:,1)==i,1)) && scn.cellID(scn.cellID(:,1)==i,6)<.05
                respPlCl(i) = 1;
            end
        end
        
        % read out special responders
        a = NaN(size(resp));
        for i = 1:size(respSpec,1)
            a(i,respSpec(i,:)==1,:) = resp(i,respSpec(i,:)==1,:);
        end
        respSpec = a; clear a;
        
        % restrict to responding cells
        respCanr = respCa(cellIDnr~=0,:,:);
        respnr = resp(cellIDnr~=0,:,:);
        respPlCl = resp(respPlCl==1 & cellIDnr~=0,:,:);
        maxAmppostnr = maxAmppost(cellIDnr~=0,:);
        maxAmpprenr = maxAmppre(cellIDnr~=0,:);
        numrespnr = numresp(cellIDnr~=0,:);
        cellIDnr = cellIDnr(cellIDnr~=0);
                
        respCa = respCa(cellID~=0,:,:);
        resp = resp(cellID~=0,:,:);
        maxAmppost = maxAmppost(cellID~=0,:);
        maxAmppre = maxAmppre(cellID~=0,:);        
        numresp = numresp(cellID~=0,:);
        cellID = cellID(cellID~=0);
        
        % cellID    firing frequency    firingfrequency during running  number of associated networks
        cellID = [cellID numresp caim.fireprop.fire(cellID,1) caim.fireprop.fire(cellID,2) caim.network.netprob(cellID)'];
        cellIDnr = [cellIDnr numrespnr caim.fireprop.fire(cellIDnr,1) caim.fireprop.fire(cellIDnr,2) caim.network.netprob(cellIDnr)'];
        
        
        
        %% Network response analysis 
        
        netresp = zeros(numstim,length(wind));
        netpost = zeros(1,numstim);
        netpre = zeros(1,numstim);
        netmean = zeros(1,length(wind));
        sumnet = zeros(numstim,length(wind));
        
        for j = 1:numstim              % Go through all onsets of stimuli to get the responses
            netresp(j,:) = caim.network.netev(find(scn.tsscn == timesscn(2*j-1))+wind);
            netmean = netmean + reshape(caim.network.netev(find(scn.tsscn == timesscn(2*j-1))+wind),size(netmean));
            netpost(j) = sum(netresp(j,find(wind==0,1):end));
            netpre(j) = sum(netresp(j,1:find(wind==0,1)));
            sumnet(j,:) = caim.network.cellsum(find(scn.tsscn == timesscn(2*j-1))+wind);
        end        
        
        
        
        %% Output
        Resp.cellID = cellID;        
        Resp.resp = resp;
        Resp.respCa = respCa;
        Resp.respPlCl = respPlCl;
        Resp.respSpec = respSpec;
        Resp.cellplc = cellplc;
        Resp.maxresppost = maxAmppost;
        Resp.maxresppre = maxAmppre;
        
        Resp.nr.cellID = cellIDnr;        
        Resp.nr.resp = respnr;
        Resp.nr.respCa = respCanr;
        Resp.nr.maxresppost = maxAmppostnr;
        Resp.nr.maxresppre = maxAmpprenr;
        
        Resp.netresp = netresp;
        Resp.netmean = netmean;
        Resp.netpost = netpost;
        Resp.netpre = netpre;    
        Resp.sumnet = sumnet;
    end
    
    %%
    if isfield(caim,'bulk')
        bulkresp = zeros(numstim,length(wind));
        bulktime = zeros(numstim,length(wind));
        bulkpost = zeros(1,numstim);
        bulkpre = zeros(1,numstim);
        bulkmean = zeros(1,length(wind));
        % Go through all onsets of stimuli to get the responses
        for j = 1:numstim              
            whichtrace = 1; %1 = raw, 2 = lowpass, 3 = bandpass, 4 = highpass
            bulkresp(j,:) = caim.bulk.trace(find(scn.tsscn == timesscn(2*j-1))+wind,whichtrace);
            bulktime(j,:) = scn.tsscn(find(scn.tsscn == timesscn(2*j-1))+wind)-timesscn(2*j-1); 
            bulkmean = bulkmean + reshape(caim.bulk.trace(find(scn.tsscn == timesscn(2*j-1))+wind,whichtrace),size(bulkmean));
            bulkpost(j) = sum(bulkresp(j,find(wind==0,1):end));
            bulkpre(j) = sum(bulkresp(j,1:find(wind==0,1)));
        end
        bulkmean = 2*bulkmean/length(timesscn);
        
        Resp.bulkresp = bulkresp;
        Resp.bulktime = bulktime;
        Resp.bulkpost = bulkpost;
        Resp.bulkpre = bulkpre;
        Resp.bulkmean = bulkmean;
        % Delay between bulk peak end GC onset 
        if isfield(caim,'C') &&  ~isempty(resp)
            a = mean(sum(resp,1),2);
            a(bulktime(1,:)<-1000 | bulktime(1,:)>1000) = 0 ;
            b = diff(bulkmean,2);
            b(end:end+2) = b(end);
            b(bulktime(1,:)<-1000 | bulktime(1,:)>1000) = 0 ;
            inout = [max(bulkmean) round(abs(bulktime(1,find(b == max(b),1)) - bulktime(1,find(a == max(a),1))))];
        else
            inout = [];
        end
        Resp.bulkinout = inout;
    end
    
    %% Behaviour readout
    speed = zeros(numstim,length(wind));
    running = zeros(numstim,length(wind));
    runtime = NaN(numstim,2);
    waittime = NaN(numstim,1);
    pupil = zeros(numstim,length(wind));
    times = zeros(numstim,length(wind));        % timestamps of traces 
    for j = 1:numstim               % Go through all onsets of stimuli to get the responses
        speed(j,:) = scn.speed(find(scn.tsscn == timesscn(2*j-1))+wind);
        running(j,1:end-1) = diff(scn.running(find(scn.tsscn == timesscn(2*j-1))+wind));
        times(j,:) = scn.tsscn(find(scn.tsscn == timesscn(2*j-1))+wind);
        % if running is initiated, read out running time
        if find(running(j,wind>0)==1,1) 
            runtime(j,1) = scn.tsscn(find(scn.tsscn>timesscn(2*j-1) & scn.running==1,1));
            runtime(j,2) = scn.totdist(scn.tsscn==runtime(j,1));
            endpoint = find(scn.tsscn>runtime(j,1) & scn.running==0,1);
            if isempty(endpoint);endpoint = length(scn.tsscn);end
            runtime(j,1) = scn.tsscn(endpoint)-runtime(j,1);
            runtime(j,2) = scn.totdist(endpoint)-runtime(j,2);
            runtime(j,1) = runtime(j,1)/1000;
        end
        % if the animal is standing, read out how long
        if scn.running(scn.tsscn==timesscn(2*j-1))==0 && ~isempty(find(scn.running(scn.tsscn<timesscn(2*j-1))==1,1)) && ~isempty(min(abs(scn.tsscn(diff(scn.running(scn.tsscn<timesscn(2*j-1)))==-1)-timesscn(2*j-1))))
            waittime(j) = min(abs(scn.tsscn(diff(scn.running(scn.tsscn<timesscn(2*j-1)))==-1)-timesscn(2*j-1)))/1000;           
        end
        times(j,:) = times(j,:)-timesscn(2*j-1);
    end
    
    
    if isfield(scn,'pupil')
        pupiltemp = zscore(scn.pupil);
    else
        pupiltemp = zeros(length(scn.tsscn),1);
    end 
    for j = 1:numstim               % Go through all onsets of stimuli to get the responses
        pupil(j,:) = pupiltemp(find(scn.tsscn == timesscn(2*j-1))+wind)-pupiltemp(find(scn.tsscn == timesscn(2*j-1))+wind(1));
%             pupil{j} = zscore(scn.pupil(find(scn.tsscn == timesscn(2*j-1))+wind));
    end
    
    Resp.times = times;    
    Resp.speed = speed;
    Resp.running = running;
    Resp.runtime = runtime;
    Resp.waittime = waittime;
    Resp.pupil = pupil;
    
    
    %% reponse quantification
    
    effects = zeros(numstim,4,3); 
    for j = 1:numstim
        for i = 1:3            
%             tempwind = find(wind==0,1)+respwind+(i-2)*respwind(end);
%             tempwind = find(wind==0,1)+(0:round(.5*Fs))+(i-2)*round(.5*Fs);
            tempwind = find(wind==0,1)+(0:round(Fs))+(i-2)*round(Fs);
            if isfield(caim,'C')
%                 Sums over cells and averages over time
                effects(j,1,i) = mean(sum(resp(:,j,tempwind)));
                effects(j,2,i) = mean(sum(respnr(:,j,tempwind)));
            end
            if isfield(caim,'bulk')
                if i == 1
                    effects(j,3,i) = bulkresp(j,find(wind==0,1));
                elseif i == 2
                    effects(j,3,i) = max(bulkresp(j,tempwind));
                end
            end
            if isfield(scn,'pupil')
                effects(j,4,i) = mean(pupil(j,tempwind));
            end
        end
    end
    Resp.effects = effects;
    Resp = downSampResp(Resp);
end

function stim = stimvect(stim)
    %% Stimulus response preference vector
    if ~isfield(stim,'resp')
        return
    end
    numbin =size(stim.stimplc,2);
    cellplc = stim.cellplc(stim.cellID(:,1),:);
%     cellplc = stim.cellplc(scn.cclust(:,cclustID.airpp)<.05,:);
    cellplchist = zeros(size(cellplc,1),numbin);
    cellcm = zeros(size(cellplc,1),1);
    for i = 1:size(cellplc,1)
        cellplchist(i,:) = histcounts(cellplc(i,cellplc(i,:)>0),1:numbin+1);
        cellcm(i) = mean(cellplc(i,cellplc(i,:)>0));
    end
    [~,b] = sort(cellcm);
    cellplchist = cellplchist(b,:);
    
    thresh = 0;
    stimspc = zeros(sum(sum(cellplchist(:,1:end),2)>thresh),numbin);
    stimprob = zeros(size(stimspc));
    span = 3;
    normspc = [stim.stimplc(end-span+1:end) stim.stimplc(1:end) stim.stimplc(1:span)];
    normspc = smooth(normspc,span)*span;
    j = 0;
    for i = 1:size(cellplc,1)
        if sum(cellplchist(i,1:end))>thresh
            j = j+1;
            stimspc(j,:) = cellplchist(i,1:end);
            stimtemp = [stimspc(j,end-span+1:end) stimspc(j,:) stimspc(j,1:span)];
    %         stimtemp = stimtemp./normspc;
            stimtemp = smooth(stimtemp,span)*span;
            stimtemp(stimtemp<=1) = 0;
            stimtemp = stimtemp./normspc;
            stimtemp(normspc<=1) = 0;
            stimprob(j,:) = stimtemp(span+1:end-span);
        end
    end
    
    stimpol = cell(size(stimspc,1),3);
    Alpha = linspace(0,2*pi,numbin);
    % theta = linspace(0,2*pi,numbin);
    rho = ones(1,length(Alpha));%(1:length(theta))/length(theta);

    for i = 1:size(stimspc,1)   % loop over cells                    
        c = [];
        d = [];
        dd = [];
        b = find(stimspc(i,:));
        for jj = 1:length(b)
            c = [c; Alpha(b(jj)),rho(b(jj)), b(jj)];  
            d = [d;exp(Alpha(b(jj))*1j)*stimspc(i,b(jj))./stim.stimplc(b(jj))];      
    %         d = [d;exp(Alpha(b(jj))*1j)];
            dd = [dd; stim.stimplc(b(jj))];
        end  

    %     stimpol{i,1} = sum(dd)*sum(d)/length(d)^2; % place coding vector
        stimpol{i,1} = sum(d)/length(d); % place coding vector
        stimpol{i,2} = abs(stimpol{i,1}); % place coding vector length
        stimpol{i,3} = c; % alpha % rho % bin 
    end
    
    stim.stimspc = stimspc;
    stim.stimprob = stimprob;
    stim.stimpol = stimpol; 
end

function stimplot(stim)
    if isempty(stim.timepnt)
        return
    end
    
    figure
    hold on
    
    plot(stim.times(1,:),20*(smooth(mean(stim.speed,1))),'c','linewidth',2)
    if isfield(stim,'pupil')
        plot(stim.times(1,:),mean(stim.pupil,1),'y','linewidth',2)
    end
    if isfield(stim,'sumnet')
        plot(stim.times(1,:),mean(stim.sumnet,1),'g','linewidth',2)
        plot(stim.times(1,:),reshape(sum(sum(stim.resp,1),2),[1 size(stim.resp,3)])./size(stim.times,1)-1,'g','linewidth',2)
        plot(stim.times(1,:),reshape(sum(sum(stim.nr.resp,1),2),[1 size(stim.nr.resp,3)])./size(stim.times,1)-1,'color',[0 .5 0],'linewidth',2)
    end
    if isfield(stim,'bulktime')
        plot(stim.bulktime(1,:),stim.bulkmean,'r','linewidth',2)
    end
    
    % plot(stim.times{1}(1:end-1),sum(cell2mat(stim.running),2)/10,'c','linewidth',2)
    plot([0 0],[-1 1],'b')
    hold off
    grid on
    axis tight
    pause(.1)
end

function Resp = downSampResp(Resp)
    Fs = 1000/scn.tsscn(end)*length(scn.tsscn); 
    FsN = 15.02;
    if Fs < 14
        newtimes = zeros(size(Resp.times,1),size(Resp.times,2)*2-1);
        for i = 1:size(Resp.times,1)
            newtimes(i,:) = Resp.times(i,1):.5*1000*(1/Fs):Resp.times(i,end); 
        end
        if showres;disp('Stimulus responses are sampled up');end
    elseif Fs > 16
        newtimes = zeros(size(Resp.times,1),ceil(size(Resp.times,2)*FsN/Fs));
        for i = 1:size(Resp.times,1)
            newtimes(i,:) = Resp.times(i,1):1000*(1/FsN):Resp.times(i,end); 
        end
        if showres;disp('Stimulus responses are sampled down');end
    else
        return
    end
    %%
     
    speed = zeros(size(newtimes));
    running = zeros(size(newtimes,1),size(newtimes,2)-1);
    pupil = zeros(size(newtimes));

    for i = 1:size(Resp.times,1)
        [~,index] = unique(Resp.times(i,:));
        speed(i,:) = interp1(Resp.times(i,index),Resp.speed(i,index),newtimes(i,:),'spline');
        running(i,:) = interp1(Resp.times(i,index),Resp.running(i,index),newtimes(i,1:end-1),'linear');
        pupil(i,:) = interp1(Resp.times(i,index),Resp.pupil(i,index),newtimes(i,:),'spline');
    end
    
    
    
    
    if isfield(caim,'C')
        
        resp = zeros(size(Resp.resp,1),size(Resp.resp,2),size(newtimes,2));
        respCa = zeros(size(Resp.resp,1),size(Resp.resp,2),size(newtimes,2));
        respnr = zeros(size(Resp.nr.resp,1),size(Resp.nr.resp,2),size(newtimes,2));
        respCanr = zeros(size(Resp.nr.resp,1),size(Resp.nr.resp,2),size(newtimes,2));
        respPlCl = zeros(size(Resp.respPlCl,1),size(Resp.respPlCl,2),size(newtimes,2));
        respSpec = zeros(size(Resp.respSpec,1),size(Resp.respSpec,2),size(newtimes,2));
        netresp = zeros(size(Resp.resp,2),size(newtimes,2));
        
        sumnet = zeros(size(Resp.resp,2),size(newtimes,2));
        
        %%
        for i = 1:size(Resp.resp,2)
            for j = 1:size(Resp.resp,1)
                resp(j,i,:) = interp1(Resp.times(i,:),reshape(Resp.resp(j,i,:),size(Resp.times(i,:))),newtimes(i,:),'next');
                respCa(j,i,:) = interp1(Resp.times(i,:),reshape(Resp.respCa(j,i,:),size(Resp.times(i,:))),newtimes(i,:),'spline');
            end
            for j = 1:size(Resp.nr.resp,1)
                respnr(j,i,:) = interp1(Resp.times(i,:),reshape(Resp.nr.resp(j,i,:),size(Resp.times(i,:))),newtimes(i,:),'next');
                respCanr(j,i,:) = interp1(Resp.times(i,:),reshape(Resp.nr.respCa(j,i,:),size(Resp.times(i,:))),newtimes(i,:),'spline');
            end
            for j = 1:size(Resp.respPlCl,1)
                respPlCl(j,i,:) = interp1(Resp.times(i,:),reshape(Resp.respPlCl(j,i,:),size(Resp.times(i,:))),newtimes(i,:),'next');
            end
            for j = 1:size(Resp.respSpec,1)
                respPlCl(j,i,:) = interp1(Resp.times(i,:),reshape(Resp.respSpec(j,i,:),size(Resp.times(i,:))),newtimes(i,:),'next');
            end
            netresp = interp1(Resp.times(i,:),Resp.netresp(i,:),newtimes(i,:),'linear');
            sumnet = interp1(Resp.times(i,:),Resp.sumnet(i,:),newtimes(i,:),'linear');         
        end
        
        netmean = interp1(Resp.times(1,:),Resp.netmean(1,:),newtimes(1,:),'linear');
        
    end
    if isfield(caim,'bulk') 
        bulkresp = zeros(size(Resp.bulkresp,2),size(newtimes,2));
        bulktime = newtimes;
%         bulkpost;
%         bulkpre;        
        %%
        for i = 1:size(Resp.times,1)
            [~,index] = unique(Resp.times(i,:));
            bulkresp = interp1(Resp.times(i,index),Resp.bulkresp(i,index),newtimes(i,:),'spline');                    
        end
        [~,index] = unique(Resp.times(1,:));
        bulkmean = interp1(Resp.times(1,index),Resp.bulkmean(1,index),newtimes(1,:),'linear');
        
        %%
        Resp.bulkresp = bulkresp(:,1:end-2);
        Resp.bulktime = bulktime(:,1:end-2);
%         Resp.bulkpost = bulkpost;
%         Resp.bulkpre = bulkpre;
        Resp.bulkmean = bulkmean(:,1:end-2);
    end


    %%
    if Fs < 14
        sht = 2;
        Resp.times = newtimes(:,1:end-sht);
        Resp.speed = speed(:,1:end-2);
        Resp.running = running(:,1:end-sht);
        Resp.pupil = pupil(:,1:end-sht);
        if isfield(caim,'C')
            
            Resp.resp = resp(:,:,1:end-sht);
            Resp.respCa = respCa(:,:,1:end-sht);
            Resp.nr.resp = respnr(:,:,1:end-sht);
            Resp.nr.respCa = respCanr(:,:,1:end-sht);
            Resp.respPlCl = respPlCl(:,:,1:end-sht);
            Resp.respSpec = respSpec(:,:,1:end-sht);
            Resp.netresp = netresp(:,1:end-sht);
            Resp.netmean = netmean(:,1:end-sht);
        %         Resp.netpost = netpost;
        %         Resp.netpre = netpre;    
            Resp.sumnet = sumnet(:,1:end-sht);
        end
        if isfield(caim,'bulk') 
            Resp.bulkresp = bulkresp(:,1:end-sht);
            Resp.bulktime = bulktime(:,1:end-sht);
        %         Resp.bulkpost = bulkpost;
        %         Resp.bulkpre = bulkpre;
            Resp.bulkmean = bulkmean(:,1:end-sht);
        end
    end
    
    if Fs > 16
        enl = 0;
        Resp.times = [newtimes(:,1:end) ones(1,enl).*newtimes(:,end)];
        Resp.speed = [speed(:,1:end) ones(1,enl).*speed(:,end)];
        Resp.running = [running(:,1:end) ones(1,enl).*speed(:,end)];
        Resp.pupil = [pupil(:,1:end) ones(1,enl).*pupil(:,end)];
        if isfield(caim,'C')
            Resp.resp = cat(3,resp(:,:,1:end), ones(1,1,enl).*resp(:,:,end));
            Resp.respCa = cat(3,respCa(:,:,1:end), ones(1,1,enl).*respCa(:,:,end));        
            Resp.nr.resp = cat(3,respnr(:,:,1:end), ones(1,1,enl).*respnr(:,:,end));
            Resp.nr.respCa = cat(3,respCanr(:,:,1:end), ones(1,1,enl).*respCanr(:,:,end));
            Resp.respPlCl = cat(3,respPlCl(:,:,1:end), ones(1,1,enl).*respPlCl(:,:,end));
            Resp.respSpec = cat(3,respSpec(:,:,1:end), ones(1,1,enl).*respSpec(:,:,end));

            Resp.netresp = [netresp(:,1:end) ones(1,enl).*netresp(:,end)];
            Resp.netmean = [netmean(:,1:end) ones(1,enl).*netmean(:,end)];
        %         Resp.netpost = netpost;
        %         Resp.netpre = netpre;    
            Resp.sumnet = [sumnet(:,1:end) ones(1,enl).* sumnet(:,end)];
        end
        if isfield(caim,'bulk') 
            Resp.bulkresp = [bulkresp(:,1:end) ones(1,enl).*bulkresp(:,end)];
            Resp.bulktime = [bulktime(:,1:end) ones(1,enl).*bulktime(:,end)];
        %         Resp.bulkpost = bulkpost;
        %         Resp.bulkpre = bulkpre;
            Resp.bulkmean = [bulkmean(:,1:end) ones(1,enl).*bulkmean(:,end)];
        end
    end
end


end


