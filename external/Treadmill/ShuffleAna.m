function [shuffle] = ShuffleAna(belt,caim,scn,numit,runshuf)
%%
if ~isfield(caim,'S')
    shuffle.shuffleDist = [];
    shuffle.pcell = [];
    shuffle.pnet = [];
    shuffle.netthresh = [];
    shuffle.networkshuffle = [];
    shuffle.shufflep = [];
    return
end

%% check if shuffled Data already exists 
if isfield(caim,'shuffle') && runshuf == 0 && numit <= size(caim.shuffle.shuffleDist,1) &&  size(caim.shuffle.cellshuffle,1) == size(caim.S_bin,1)
    disp('loaded shuffled data')
    shuffle = caim.shuffle;
    return
end
%%

disp('run shuffle analysis')
thresh = 3:10;
networkshuffle = zeros(numit,3,length(thresh));
fireshuffle = zeros(numit,3);
bulkshuffle = zeros(numit,1);
bulkstimshuffle = zeros(numit,1);
cellshuffle = zeros(size(caim.S_bin,1),5,numit);
parfor i = 1:numit   
    %%   
    caim2 = ShuffleData(caim);    
    nstemp = zeros(3,length(thresh));
    for j = 1:length(thresh)
        caim2.netthresh = thresh(j);
        network = netevents(caim2,scn,0);
        if ~isempty(network.netfreq)
            nstemp(:,j) = network.netfreq(1,:);
        end
    end
    networkshuffle(i,:,:) = nstemp;
    %%
    fireprop = firefreq(caim2,scn);      
    celltemp = zeros(size(caim2.S_bin,1),5);        
    fireshuffle(i,:) = fireprop.meanfire(2,:);   
    celltemp(:,1) = sum(caim2.S_bin(:,scn.running==1),2)./sum(caim2.S_bin,2);

    if isfield(scn,'airpuff')
        scn1 = stimcor(belt,caim2,scn,0);
        stim = scn1.airpuff;

        % cell is a responder
        celltemp(stim.cellID(:,1),2) = 1;
        % cell is a non responder 
        celltemp(stim.nr.cellID(:,1),3) = 1;
        % number of stimulus related events       
        celltemp(stim.cellID(:,1),4) = stim.cellID(:,2);
        celltemp(stim.nr.cellID(:,1),4) =stim.nr.cellID(:,2);
        % number of reponded stimuli
        celltemp(stim.cellID(:,1),5) = stim.cellID(:,3);
        celltemp(stim.nr.cellID(:,1),5) = stim.nr.cellID(:,3);

    else
        celltemp(:,2:5) = NaN;
    end


    cellshuffle(:,:,i) = celltemp;

    %%
    if isfield(caim,'bulk')
        belt1 = ShuffleData(belt);
        caim2 = ShuffleData(caim,1);
        scn1 = stimcor(belt1,caim2,scn,0);

        win = find(abs(scn.network.times(1,:))==min(abs(scn.network.times(1,:)))):find(abs(scn.network.times(1,:))==min(abs(scn.network.times(1,:))))+15;
        a = max(scn1.network.bulkmean(win));
        bulkshuffle(i) = a(1);
        if isfield(scn,'airpuff')
            bulkstimshuffle(i) = max(scn1.airpuff.bulkmean(win));
        end
    end
end   

%% fit network p

pnet = zeros(3,length(thresh));
a = zeros(length(thresh),3);
b = zeros(length(thresh),3);
for i = 1:length(thresh)

    caim.netthresh = thresh(i);
    network = netevents(caim,scn,0);
    a(i,:) = network.netfreq(1,:);
    b(i,:) = mean(networkshuffle(:,:,i));
    
    pnet(1,i) = sum(a(i,1)<=networkshuffle(:,1,i))/numit;
    pnet(2,i) = sum(a(i,2)<=networkshuffle(:,2,i))/numit;
    pnet(3,i) = sum(a(i,3)<=networkshuffle(:,3,i))/numit;
    for j = 1:3
        if i>1 && pnet(j,i)>=.05 && pnet(j,i)>pnet(j,i-1)
            pnet(j,1:i-1) = 1;
        end
    end
end

pnet = [pnet; thresh];
a = find(pnet(1,:)<.05,1)+1;
if isempty(a);a = 1;end
thresh = thresh(a);


%% comparisson values

fireprop = firefreq(caim,scn);
scn = stimcor(belt,caim,scn,0);

celltemp(:,1) = sum(caim.S_bin(:,scn.running==1),2)./sum(caim.S_bin,2);

if isfield(scn,'airpuff')
    stim = scn.airpuff;

    % number of stimulus related events
    celltemp(stim.cellID(:,1),4) = stim.cellID(:,2);
    celltemp(stim.nr.cellID(:,1),4) =stim.nr.cellID(:,2);

    % number of responded stimuli
    celltemp(stim.cellID(:,1),5) = stim.cellID(:,3);
    celltemp(stim.nr.cellID(:,1),5) = stim.nr.cellID(:,3);
else
    celltemp(:,2:5) = NaN;
end

if isfield(caim,'bulk')  
    win = find(min(abs(scn.network.times(1,:)))==abs(scn.network.times(1,:))):find(min(abs(scn.network.times(1,:)))==abs(scn.network.times(1,:)))+15; 
    bulkmean = max(scn.network.bulkmean(win));
    if isfield(scn,'airpuff')
        bulkstimmean = max(scn.airpuff.bulkmean(win));
    end
end




%%
p = pnet(:,a);
p(4) = sum(fireshuffle(:,2)>fireprop.meanfire(2,2))/numit;
p(5) = sum(fireshuffle(:,3)>fireprop.meanfire(2,3))/numit;
if isfield(caim,'bulk')  
    p(6) = sum(bulkshuffle>bulkmean)/numit;
    if isfield(scn,'airpuff')
        p(7) = sum(bulkstimshuffle>bulkstimmean)/numit;
    else
        p(7) = NaN;
    end
else
    p(6) = NaN;
    p(7) = NaN;
end

shuffleDist = [networkshuffle(:,:,a) fireshuffle bulkshuffle bulkstimshuffle];

%% individual cell p values
pcell = zeros(size(celltemp));
for i = 1:size(pcell,1)
    % fraction of running related events
    pcell(i,1) = sum(cellshuffle(i,1,:)>celltemp(i,1))/numit;
    % fraction of cell being a random responder
    pcell(i,2) = sum(cellshuffle(i,2,:))/numit;
    % fraction of cell being a random non responder
    pcell(i,3) = sum(cellshuffle(i,3,:))/numit;
    % probability of having more response events by chance
    pcell(i,4) = sum(cellshuffle(i,4,:)>celltemp(i,4))/numit;
    % probability of having more responded stimuli by chance
    pcell(i,5) = sum(cellshuffle(i,5,:)>celltemp(i,5))/numit;
end

%% outpout

shuffle.shuffleDist = shuffleDist;
shuffle.cellshuffle = cellshuffle;
shuffle.pcell = pcell;
shuffle.pnet = pnet;
shuffle.netthresh = thresh;
shuffle.networkshuffle = networkshuffle;
shuffle.shufflep = p;



%%
% figure
% subplot(3,3,1)
% histogram(networkshuffle(:,1),'Normalization','probability')
% hold on 
% plot([network.netfreq(1,1) network.netfreq(1,1) ],[0 .4])
% 
% subplot(3,3,2)
% histogram(networkshuffle(:,2),'Normalization','probability')
% hold on 
% plot([network.netfreq(1,2) network.netfreq(1,2) ],[0 .4])
% 
% subplot(3,3,3)
% histogram(networkshuffle(:,3),'Normalization','probability')
% hold on 
% plot([network.netfreq(1,3) network.netfreq(1,3) ],[0 .4])
% 
% subplot(3,3,5)
% histogram(fireshuffle(:,2),'Normalization','probability')
% hold on
% plot([fireprop.meanfire(2,2) fireprop.meanfire(2,2)],[0 0.2])
% 
% subplot(3,3,6)
% histogram(fireshuffle(:,3),'Normalization','probability')
% hold on
% plot([fireprop.meanfire(2,3) fireprop.meanfire(2,3)],[0 0.2])
% 
% subplot(3,3,7)
% histogram(bulkshuffle,'Normalization','probability')
% hold on
% plot([bulkmean bulkmean],[0 0.2])
% 
% subplot(3,3,8)
% histogram(bulkstimshuffle,'Normalization','probability')
% hold on
% plot([bulkstimmean bulkstimmean],[0 0.2])


end

function DataOut = ShuffleData(DataIn,netshuf)

if isfield(DataIn,'S') && nargin < 2
    caim = DataIn;
    S_bin = caim.S_bin;
    % shuffle everything
    for j = 1:size(S_bin,1)
        a = randperm(size(S_bin,2));
        caim.S_bin(j,:) = S_bin(j,a);
        caim.C(j,:) = caim.C(j,a);
    end
       
    % shift traces randomly with respect to each other
    %     for j = 1:size(S_bin,1)
    %         a = randperm(size(S_bin,2),1);
    %         S_temp(j,:,1) = S_bin(j,[a:end 1:a-1]);
    %         S_temp(j,:,2)  = C(j,[a:end 1:a-1]);
    %     end
    
%     caim.shuffle = 1;
    DataOut = caim;
    if isfield(caim,'shuffle')
        caim  = rmfield(caim,'shuffle');
    end
elseif isfield(DataIn,'S') && netshuf == 1
    caim = DataIn;
    a = randperm(size(caim.network.netev,1));
    caim.network.netev = caim.network.netev(a);
%     caim.shuffle = 1; 
    DataOut = caim;
    if isfield(caim,'shuffle')
        caim  = rmfield(caim,'shuffle');
    end
elseif isfield(DataIn,'time')
    belt = DataIn;
    a = randperm(size(belt.airpuff,1));
    belt.airpuff = belt.airpuff(a);
    DataOut = belt;    
end
end
