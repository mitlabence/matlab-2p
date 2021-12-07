function cellraster(CAIMcorr,mouse,expID)
%%
load('cclustID.mat','cclustID');
meanf = cclustID.meanf;             % mean frequency of cell firing
plcfld =  cclustID.plcfld;          % placefield center
% plcfldp = cclustID.plcfldp;       % placefield p value
plclength = cclustID.plclength;     % size of place field 
plcvang = cclustID.plcvctang;       % place coding vector angle
plcvct = cclustID.plcvct;           % place coding vector length
plcvctp = cclustID.plcvctp;         % place coding vector p value
nstim = cclustID.nstim;             % # responded stimuli
netperc = cclustID.netperc;         % percentage of network related firing
netprob = cclustID.netprob;         % netprob
mpp = cclustID.mpp;                 % mean MPP input
runperc = cclustID.runperc;         % percentage of running associated cells
runnerp = cclustID.runnerp;         % p value running
airpp = cclustID.airpp;             % p value number of response events
followcell = cclustID.followcell;   % logical if cell was in FOV of all experiments

%%
emptyCAIM = [];
for i = 1:length(expID)
    if isempty(CAIMcorr(expID(i),mouse).A)
        emptyCAIM = [emptyCAIM i];
    end
end
expID(emptyCAIM) = [];

netraster = cell(length(expID),1);
netnum = zeros(length(expID)+1,1);
isplace = cell(length(expID),1);
isstim = cell(length(expID),1);
isboth = cell(length(expID),1);
isthere = cell(length(expID),1);

for i = 1:length(expID)
    trial = expID(i);
    netraster{i} = CAIMcorr(trial,mouse).network.netraster;
    netnum(i+1) = size(netraster{i},2)+netnum(i);
    cclusttemp = CAIMcorr(trial,mouse).cclust;
%     active field sorting
%     netraster{i} = netraster{i}(CAIMcorr(trial,mouse).inField==1,:);
%     cclusttemp = cclusttemp(CAIMcorr(trial,mouse).inField==1,:);
    netraster{i}(CAIMcorr(trial,mouse).inField==0,:) = 0;
    cclusttemp(CAIMcorr(trial,mouse).inField==0,:) = 0;

    isplace{i} = (cclusttemp(:,plcvct)>0 & cclusttemp(:,plcvctp)<=.05) & ~(cclusttemp(:,nstim)>4 & cclusttemp(:,airpp)<=.05);
    isstim{i}  = ~(cclusttemp(:,plcvct)>0 & cclusttemp(:,plcvctp)<=.05) & (cclusttemp(:,nstim)>4 & cclusttemp(:,airpp)<=.05);
    isboth{i}  = (cclusttemp(:,plcvct)>0 & cclusttemp(:,plcvctp)<=.05) & (cclusttemp(:,nstim)>4 & cclusttemp(:,airpp)<=.05);
    isthere{i} = ~isnan(cclusttemp(:,1));
end

% fill up individual indicators so that they all have smae length
for i = 1:length(expID-1)
    netraster{i}(size(netraster{end},1),:) = 0;
    isplace{i}(size(netraster{end},1)) = 0;
    isstim{i}(size(netraster{end},1)) = 0;
    isboth{i}(size(netraster{end},1)) = 0;
    isthere{i}(size(netraster{end},1)) = 0;
end

netraster = cell2mat(netraster');

% transfer the individual cells to one matrix with aditional indicator in
% first coloumn
isplace = cell2mat(isplace');
isplace(:,2:end+1) = isplace; 
isplace(:,1) = sum(isplace(:,2:end),2);
isplace(isplace>1) = 1;

isstim = cell2mat(isstim');
isstim(:,2:end+1) = isstim; 
isstim(:,1) = sum(isstim(:,2:end),2);
isstim(isstim>1) = 1;

isboth = cell2mat(isboth');
isboth(:,2:end+1) = isboth; 
isboth(:,1) = sum(isboth(:,2:end),2);
isboth(isboth>1) = 1;

isstim(isboth(:,1)==1,1)=0;
% isstim(isplace(:,1)==1,1)=0;
isplace(isboth(:,1)==1,1)=0;
isplace(isstim(:,1)==1,1)=0;

isthere = cell2mat(isthere');
isthere(:,2:end+1) = isthere;

newsort = [find(isplace(:,1)); find(isboth(:,1)); find(isstim(:,1))];
isplace = isplace(newsort,:);
isstim = isstim(newsort,:);
isboth = isboth(newsort,:);
isthere = isthere(newsort,:);
netraster = netraster(newsort,:);


% figure('color',[1 1 1],'position',[500 50 1.5*[420 594]],'visible','on')
hold on

for i = 1:length(netnum)-2
    plot([netnum(i+1) netnum(i+1)],[0 size(isplace,1)] ,':','color',[0 0 0])
end

if ~isempty(find(isboth,1))
    plot([0 size(netraster,2)],[find(isboth(:,1),1)-.5 find(isboth(:,1),1)-.5] ,':','color',[0 0 0])
end

if ~isempty(find(isstim,1))
    plot([0 size(netraster,2)],[find(isstim(:,1),1)-.5 find(isstim(:,1),1)-.5] ,':','color',[0 0 0])
end

sumplace = zeros(size(netraster,2),1);
sumstim = zeros(size(netraster,2),1);
sumboth = zeros(size(netraster,2),1);
sumnone = zeros(size(netraster,2),1);
for j = 2:length(netnum)
    win = netnum(j-1)+1:netnum(j); 
    % lifelines 
    actcells = find(isthere(:,j));    
    for i = 1:length(actcells)
        plot([win(1) win(end)],[actcells(i) actcells(i)],'color',[.8 .8 .8])
    end
    % placecells
    evcells = find(isplace(:,j));
    for i = 1:length(evcells)
        events = find(netraster(evcells(i),win)); 
        scatter(events+netnum(j-1),ones(1,length(events))*evcells(i),'filled','MarkerFaceColor',[0 176/255,240/255]);
        sumplace(win) = sumplace(win) + netraster(evcells(i),win)';
    end
    % AP cells
    evcells = find(isstim(:,j));
    for i = 1:length(evcells)
        events = find(netraster(evcells(i),win)); 
        scatter(events+netnum(j-1),ones(1,length(events))*evcells(i),'filled','MarkerFaceColor',[1 0 1]);  
        sumstim(win) = sumstim(win) + netraster(evcells(i),win)';
    end
    % both 
    evcells = find(isboth(:,j));
    for i = 1:length(evcells)
        events = find(netraster(evcells(i),win)); 
        scatter(events+netnum(j-1),ones(1,length(events))*evcells(i),'filled','MarkerFaceColor',[1 1 0]);  
        sumboth(win) = sumboth(win) + netraster(evcells(i),win)';
    end
    % none
    evcells = find(~isstim(:,j)&~isplace(:,j)&~isboth(:,j));
    for i = 1:length(evcells)
        events = find(netraster(evcells(i),win)); 
        scatter(events+netnum(j-1),ones(1,length(events))*evcells(i),'filled','MarkerFaceColor',[.5 .5 .5]);  
        sumnone(win) = sumnone(win) + netraster(evcells(i),win)';
    end    
end

sumnone(:) = 0;
p = bar(-[sumstim sumboth sumplace sumnone],'stacked');
piecol = {[1 0 1],[1 1 0],[0 176/255,240/255],[.5 .5 .5]};

for i = 1:4
    pp = p(i);
    pp.FaceColor  = piecol{i};
end
axis tight

end
