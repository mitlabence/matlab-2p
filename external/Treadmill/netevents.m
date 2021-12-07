function network = netevents(caim,scn,showres) 
%%
if isstruct(caim) && isfield(caim,'S_bin')
    S_bin = caim.S_bin;
    C = caim.C;  
else
    network = [];
    return
end

if nargin < 3
    showres = 1;
end

%% Detect Network events 

Fs = 1000*length(scn.tsscn)/scn.tsscn(end);     % Sampling frequency
T = 1/Fs;                                       % Sampling period
win = .2;                                       % binning window in s
wind = 1;                                       % window in between two events to be considered as one
span = round(win/T);                            % Span of moving average window
if isfield(caim,'netthresh')
    thresh = caim.netthresh;
else
    thresh = ceil(0.01 * size(S_bin,1));        % minimal number of cells needed for event
    if thresh < 4;thresh = 4;end                % lower bound for threshold (original 4)
end

netev = smooth(sum(S_bin),span)*span;
% netev = sum(S_bin);
netval = netev;
netev(netev<=thresh) = 0;
netev(netev>thresh) = 1;
wind = round(wind/T);
netev(end-span:end) = 0;
netev(1:1+span) = 0;
%%
for j = 1:2
    %%
    
    stepup = find(diff(netev)==1);
    stepdown = find(diff(netev)==-1); 
    if length(stepup) < length(stepdown)
        stepdown = stepdown(2:end);
    end
    netpos = zeros(1,length(stepup));
    for i = 1 : length(stepdown)-1
        if stepup(i+1)-stepdown(i) < wind
            netev(stepdown(i):stepup(i+1)) = 1;
        end
        netpos(i) = stepup(i)+find(netval(stepup(i):stepdown(i)) == max(netval(stepup(i):stepdown(i))),1);
    end   
end


netpos(i+1) = stepup(i+1)+find(netval(stepup(i+1):stepdown(i+1)) == max(netval(stepup(i+1):stepdown(i+1))),1);
if isempty(netpos) || netpos(1) == 0  
    if showres
        disp('No Network events detected')
    end
    network.thresh = [];
    network.netpos = [];
    network.netfreq = zeros(4,3);
    network.netev = netev;
    network.netID = [];
    network.netprob = zeros(1,size(caim.Y,1));
    network.netnum = [];
    network.replay = [];
    network.coact = [];
    network.netresp = zeros(5,150);
    network.cellsum = zeros(size(netev));
    network.intEvInt = zeros(1,200);
    network.netfirecorr = [];
    return
end
netpos = netpos-1;

%% Cell ID per event

netID = cell(1,length(netpos));
netprob = zeros(1,size(S_bin,1));
netnum = zeros(1,length(netpos));
for i = 1 : length(netpos)
    a = zeros(1,size(S_bin,1));
    for j = 1:size(S_bin,1)
        if find(S_bin(j,netpos(i)-span:netpos(i)+span),1)
            a(j) = 1;
        end
    end
    if length(find(a)) > thresh
        netID{i} = find(a);
        netnum(i) = length(find(a));
        netprob = netprob +a;  
    end
end

netpos = netpos(netnum>0);
netID = netID(netnum>0);
netnum = netnum(netnum>0);

netev(:) = 0;
netev(netpos) = 1;

%% calculate frequencies in different stages

lges = scn.tsscn(end)/1000/60;                     %total time in min
lrun = T*length(scn.running(scn.running==1))/60;   %time running
lstand = T*length(scn.running(scn.running==0))/60; %time standing

netfreq = zeros(2,3);

netfreq(1,1) = sum(netev);
netfreq(1,2) = sum(netev(scn.running==1));
netfreq(1,3) = sum(netev(scn.running==0));

netfreq(2,1) = netfreq(1,1)/lges;
netfreq(2,2) = netfreq(1,2)/lrun;
netfreq(2,3) = netfreq(1,3)/lstand;

%% The raster
netraster = zeros(size(S_bin,1),length(netID));
for i = 1:length(netID)
    netraster(netID{i},i)=1;
end

%% filter all events that occured during running
standnet = scn.running(netpos)==0;
netID = netID(standnet);
netpos = netpos(standnet);
netnum = netnum(standnet);
netraster = netraster(:,standnet);
netprob = sum(netraster,2)';
netev(:) = 0;
netev(netpos) = 1;
%%
if showres
    disp([num2str(length(netpos)) ' Network events detected'])
end
%% inter event interval
if ~isempty(find(scn.running,1))
    start = find(diff(scn.running)==1);
    stop = find(diff(scn.running) == -1);
    start = reshape(start,1,length(start));
    stop = reshape(stop,1,length(stop));
    if stop(1)>start(1)
        stop = [1 stop];
    end
    if start(end)<stop(end)
        start = [start length(netev)];
    end

    c = [];
    for i = 1 :length(stop)
        a = netev(stop(i):start(i));
        time = scn.tsscn(stop(i):start(i));
        b = time(a == 1);
        if length(b) > 2
            c = [c; diff(b)];
        end
    end
    intEvInt = histcounts(c,0:100:20000);
else
    intEvInt = NaN(1,200);
end
%%

cellsum = sum(C);
cellsum = debleach(scn.tsscn/1000,cellsum );
% cellsum = cellsum';
cellsum = zscore(cellsum);

netfreq(3,1) = mean(cellsum);
netfreq(3,2) = mean(cellsum(scn.running==1));
netfreq(3,3) = mean(cellsum(scn.running==0));

netfreq(4,1) = std(cellsum);
netfreq(4,2) = std(cellsum(scn.running==1));
netfreq(4,3) = std(cellsum(scn.running==0));

%% replay

if isfield(scn,'cellID')
    % plccell = scn.cellID(scn.cellID(:,2)>0 & scn.cellID(:,3)<=0.05,1); % Dombeck
    plccell = scn.cellID(scn.cellID(:,6)<=0.05,1); % loso
    plchist = zeros(1,length(netID));

    for i = 1:length(netID)
        k = 0;
        for j = 1:length(plccell)
            if find(plccell(j) == netID{i})
                k = k+1;
            end
        end
        plchist(i) = k;
    end
    replay = histcounts(plchist,-.5: 1 : 9.5);
else
    replay = NaN(1,10);
end

%% Read out activity of participants and non-participants
cellnum = size(C,1);
wind = [-2 8];                                                        % size of window (in s) to look for response
wind = round(wind(1)*Fs):round(wind(end)*Fs);
if length(wind)==151;wind=wind(1:150);end
netpostemp = netpos(netpos>-wind(1)&size(C,2)-netpos>wind(end));
resp = nan(cellnum,length(netnum),length(wind));                  % binary Data to look for responses
respCa = nan(cellnum,length(netnum),length(wind));                % traces of Data to look for responses        
respnr = nan(cellnum,length(netnum),length(wind));                  % binary Data to look for responses
respnrCa = nan(cellnum,length(netnum),length(wind));                % traces of Data to look for responses        

for j = 1:length(netpostemp)
    cellID = netID{j};
    for i = 1:length(cellID)        
        resp(cellID(i),j,:) = S_bin(cellID(i),netpostemp(j)+wind);
        respCa(cellID(i),j,:) = C(cellID(i),netpostemp(j)+wind);
    end
    cellID = 1:cellnum;
    cellID(netID{j}) = [];
    for i = 1:length(cellID)        
        respnr(cellID(i),j,:) = S_bin(cellID(i),netpostemp(j)+wind);
        respnrCa(cellID(i),j,:) = C(cellID(i),netpostemp(j)+wind);
    end
end

resp = [wind/Fs;
    permute(nanmean(nanmean(resp,1),2),[1 3 2]);
    permute(nanmean(nanmean(respnr,1),2),[1 3 2]);
    permute(nanmean(nanmean(respCa,1),2),[1 3 2]);
    permute(nanmean(nanmean(respnrCa,1),2),[1 3 2])];

resp = downSampResp(resp);

%% check for co-activity between cells
coact = zeros(cellnum);
for i = 1:length(netID)
    a = netID{i};
    for j = 1:length(a)
        k = j+1;
        while k<=length(a)
            coact(a(j),a(k)) = coact(a(j),a(k))+1;
            coact(a(k),a(j)) = coact(a(k),a(j))+1;           
            k = k+1;
        end
    end
end

% coact = coact./length(netID);

%% cos-distance between network events

netdist = zeros(size(netraster,2));
for i = 1:length(netdist)
    for j = 1:length(netdist)
        netdist(i,j) = (netraster(:,j)'*netraster(:,i))./(sqrt(sum(netraster(:,j)))* sqrt(sum(netraster(:,i))) );
    end
end
%% cos-distance between cells

celldist = zeros(size(netraster,1));
for i = 1:length(celldist)
    for j = 1:length(celldist)
        celldist(i,j) = (netraster(j,:)*netraster(i,:)')./(sqrt(sum(netraster(j,:)))* sqrt(sum(netraster(i,:))) );
    end
end
%% Check out places of net events for all cells

b = find(netev);
netplace = zeros(cellnum,150);
for i = 1:size(C,1)
    a = netraster(i,:)==1;
    c = b(a);
    netplace(i,:) = histcounts(scn.distance(c),0:10:1500);
end

%% Read out distance of PCs to next net event

if isfield(scn,'cellID')
    cellID = scn.cellID;
    spcsort = scn.spcsort(cellID(:,6)<=.05,:,:);
    plcfield = scn.plcfield(cellID(:,6)<=.05,:,1);
    fieldstemp = scn.spcpol(cellID(:,6)<=.05,4);
    fields = zeros(size(plcfield));

    for i = 1:size(fields,1)
        fields(i,:) = fieldstemp{i};   
        fields(i,:) = mat2gray(fields(i,:));
        fields(i,:) = fields(i,:)>.66;
    end

    cellID = cellID(cellID(:,6)<=.05,:);
    b = find(netev);
    netplace = zeros(size(cellID,1),150);
    fieldnet = NaN(size(spcsort,3),8,size(cellID,1));
    % read out overall place prefrerence
    if ~isempty(cellID)
        cnet = cell(size(cellID,1),5);
        for i = 1:size(cellID,1)
            a = netraster(cellID(i,1),:)==1;
            cnet{i,1} = cellID(i,1);            % Cell ID
            cnet{i,2} = cellID(i,7);            % Place ifeld position
            cnet{i,3} = b(a);                   % indices
            cnet{i,4} = scn.tsscn(b(a));        % timepoints
            cnet{i,5} = scn.distance(b(a));     % positions
        end
    else 
        cnet = [];
    end
    % go through every round and check activity of every cell
    for i = 1:size(cellID,1)
        a = netraster(cellID(i,1),:)==1;
        cnettime = b(a);
        plcnet = [];
        for j = 1:size(spcsort,3)       
            % time of round
            thisround = scn.rounds == j;
            roundtime = scn.tsscn(thisround);
            rounddist = scn.distance(thisround);
            roundtotdist = scn.totdist(thisround);
            roundrun = scn.running(thisround);

            % activity in round
            cellact = caim.S_bin(cellID(i,1),thisround);
            act = find(thisround & caim.S_bin(cellID(i,1),:)'==1);
            act(:,2) = roundtime(cellact==1);
            act(:,3) = rounddist(cellact==1);
            act(:,4) = roundrun(cellact==1);
            for k = 1:size(act,1)
                % field crossing  
                if ceil(act(k,3)/10)>0
                    act(k,5) = fields(i,ceil(act(k,3)/10)); 
                    % check if cell fired during runnning within field
                    if ~isempty(cnettime) && act(k,4)==1 && act(k,5)==1
                        % time & place to next network before and after
                        if ~isempty(find(cnettime>act(k,1),1))
                            act(k,6) = cnettime(find(cnettime>act(k,1),1));
                            act(k,7) = scn.tsscn(act(k,6));
                            act(k,8) = scn.distance(act(k,6));
                        else
                            act(k,6:8) = NaN;
                        end

                        if ~isempty(find(cnettime<act(k,1),1,'last'))
                            act(k,9) = cnettime(find(cnettime<act(k,1),1,'last'));
                            act(k,10) = scn.tsscn(act(k,9));
                            act(k,11) = scn.distance(act(k,9));
                        else
                            act(k,9:11) = NaN;
                        end
                    else
                        act(k,6:11) = NaN;
                    end
                else
                    act(k,5) = 0;
                end
            end
            if ~isempty(act) && ceil(act(k,3)/10)>0
                plcnet = [plcnet; act(act(:,4)==1 & act(:,5)==1,:)];            
            end
            % Crossing time of place field
            fieldcross = abs(rounddist-cellID(i,7))==min(abs(rounddist-cellID(i,7))); 
            if ~isempty(fieldcross)
                fieldnet(j,1,i) = roundtime(find(fieldcross,1));
                fieldnet(j,2,i) = roundtotdist(find(fieldcross,1));
                if ~isempty(act) && sum(act(:,4)==1 & act(:,5)==1)>=1            
                    fieldcrossindex = find(thisround,1)+find(fieldcross,1);           
                    netbefore = cnettime(find(cnettime<fieldcrossindex,1,'last'));
                    if ~isempty(netbefore)
                        fieldnet(j,3,i) = scn.tsscn(netbefore);
                        fieldnet(j,4,i) = scn.totdist(netbefore);
                        fieldnet(j,7,i) = cellID(i,7);
                        fieldnet(j,8,i) = scn.distance(netbefore);
                    end
                    netafter = cnettime(find(cnettime>fieldcrossindex,1));
                    if ~isempty(netafter)
                        fieldnet(j,5,i) = scn.tsscn(netafter);
                        fieldnet(j,6,i) = scn.totdist(netafter);
                        fieldnet(j,7,i) = cellID(i,7);
                        fieldnet(j,8,i) = scn.distance(netafter);
                    end
                end   
            end
        end
    end
else
    cnet = [];
    fieldnet = [];
end

%% Correlation vector for network and individual firing

if isfield(caim,'fireprop') %%&& isfield(caim.network,'netprob')
    netfirecorr = [caim.fireprop.fire(netprob~=0,1)';netprob(netprob~=0)/lges];
else
    netfirecorr = [];
end

%% Mean of shuffled network events
if isfield(caim,'shuffle') % && caim.shuffle ~=1
    shuffleDist = caim.shuffle.shuffleDist;
    ShufMean = [mean(shuffleDist(:,1:3),1);std(shuffleDist(:,1:3),1)];
    pfit = caim.shuffle.pnet;
    netp = caim.shuffle.shufflep;
else
    ShufMean = NaN(1,3);  
    pfit = [];
    netp = [];
end
%%
network.thresh = thresh;
network.pfit = pfit;
network.netp = netp;
network.netev = netev;
network.netpos = scn.tsscn(netpos);
network.netID = netID;
network.netprob = netprob;
network.netfreq = netfreq;
network.netraster = netraster;
network.netdist = netdist;
network.celldist = celldist;
network.netnum = netnum;
network.netresp = resp;
network.cellsum = cellsum;
network.replay = replay;
network.coact = coact;
network.intEvInt = intEvInt;
network.netplace = netplace;
network.netfirecorr = netfirecorr;
network.fieldnet = fieldnet;
network.cnet = cnet;
network.ShufMean = ShufMean;
%% Correlation to running and fft

% d = scn.tsscn/1000;
% a = sum(C);
% a = debleach(d,a);
% a = a-min(a);
% a = a/max(a);
% speed = diff(scn.distance);
% speed(speed<-20) = 0;
% speed(end+1)=speed(end);
% speed = smooth(speed,30);
% speed = speed-min(speed);
% speed = speed/max(speed);

%% Fourier transformation of summed events
% Fs = length(scn.tsscn)/scn.tsscn(end)*1000;            % Sampling frequency
% T = 1/Fs;             % Sampling period
% a1 = a(scn.running==0);
% L = length(a1);             % Length of signal
% Y = fft(a1);
% P1 = abs(Y/L);
% P1 = P1(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% f1 = Fs*(0:(L/2))/L;
% 
% 
% a2 = a(scn.running==1);
% L = length(a2);             % Length of signal
% Y = fft(a2);
% P2 = abs(Y/L);
% P2 = P2(1:L/2+1);
% P2(2:end-1) = 2*P2(2:end-1);
% f2 = Fs*(0:(L/2))/L;
% 
% plot(f1,smooth(P1,50))
% hold on
% plot(f2,smooth(P2,50))
% hold off

%%
%     
%     [P,Q] = rat(length(P2)/length(P1));
%     ynew = smooth(resample(P1,P,Q),100);
%     ynew2 = smooth(P2,100);
%     if length(ynew) > length(ynew2)
%         ynew = ynew(1:length(ynew2));
%     elseif length(ynew) < length(ynew2)
%         ynew2 = ynew(1:length(ynew));
%     end
%     ratio(1) = min(abs(ynew2(f2>.1 & f2<4)./ynew(f2>.1 & f2<4)));
%     if  Fs > 10
%         ratio(2) = max(abs(ynew2(f2>4 & f2<15)./ynew(f2>4 & f2<15)));
%     else 
%         ratio(2) = 0;
%     end
%     plot(f2,abs(ynew2./ynew))
% %     plot(f2,ynew)
% %     hold on
% %     plot(f2,ynew2)
% %     hold off



function newresp = downSampResp(resp)   
    if Fs < 14
        newresp(1,:) = resp(1):.5*(1/Fs):resp(1,end); 
        if ~isnan(resp(2,1))
            newresp(2,:) = interp1(resp(1,:),resp(2,:),newresp(1,:),'next');
            newresp(3,:) = interp1(resp(1,:),resp(3,:),newresp(1,:),'next');
            newresp(4,:) = interp1(resp(1,:),resp(4,:),newresp(1,:),'spline');
            newresp(5,:) = interp1(resp(1,:),resp(5,:),newresp(1,:),'spline');
        else
            newresp(2:5,:) = nan(4,length(newresp));
        end
        sht = abs(length(newresp(1,:))-150);
        newresp = newresp(:,1:end-sht);
        if showres;disp('Network are sampled up');end
    elseif Fs > 16
        newresp(1,:) = resp(1):2*(1/Fs):resp(1,end); 
        if ~isnan(resp(2,1))
            newresp(2,:) = interp1(resp(1,:),resp(2,:),newresp(1,:),'next');
            newresp(3,:) = interp1(resp(1,:),resp(3,:),newresp(1,:),'next');
            newresp(4,:) = interp1(resp(1,:),resp(4,:),newresp(1,:),'spline');
            newresp(5,:) = interp1(resp(1,:),resp(5,:),newresp(1,:),'spline');
        else
            newresp(2:5,:) = nan(4,length(newresp));
        end
        
        enl = abs(length(newresp(1,:))-150);
        newresp = cat(2,newresp(:,1:end), ones(1,enl).*newresp(:,end));
        if showres;disp('Network responses are sampled down');end
    else
        newresp=resp;
    end   
end
end
