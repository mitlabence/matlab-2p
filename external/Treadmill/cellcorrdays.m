function [cellcorr,TodayCorr] = cellcorrdays(CAIMcorr,mouse,expID)

experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};

load('cclustID.mat','cclustID');

%%
emptyCAIM = false(1,length(expID));
for i = 1:length(expID)
    if isempty(CAIMcorr(expID(i),mouse).A)
        emptyCAIM(i) = true;
    end
end
expID(emptyCAIM) = [];
cclust = cell(length(expID),1);
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
    netraster{i}(CAIMcorr(trial,mouse).inField==0,:) = 0;
    cclusttemp(CAIMcorr(trial,mouse).inField==0,:) = NaN;

    isplace{i} = (cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & ~(cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
    isstim{i}  = ~(cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & (cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
    isboth{i}  = (cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & (cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
    isthere{i} = ~isnan(cclusttemp(:,1));
    cclust{i} = cclusttemp;
end

%% fill up individual indicators so that they all have same length

for i = 1:length(expID-1)
    netraster{i}(size(netraster{end},1),:) = 0;
    isplace{i}(size(netraster{end},1)) = 0;
    isstim{i}(size(netraster{end},1)) = 0;
    isboth{i}(size(netraster{end},1)) = 0;
    isthere{i}(size(netraster{end},1)) = 0;
    cclust{i}(size(netraster{end},1),1) = NaN;
end

%%
cclusttemp = cell2mat(cclust(1));
for i = 2:length(cclust)
    cclusttemp = cat(3,cclusttemp,cell2mat(cclust(i)));
end
cclust = cclusttemp;

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
% isplace(isstim(:,1)==1,1)=0;

isthere = cell2mat(isthere');
isthere(:,2:end+1) = isthere;


%% correlation thru trials following cells

int = 1:length(expID);

%% sum (mean) of all correlated cells to stim cells
cellID = isthere(:,1)>0;
cclusttemp = cclust(cellID,:,:);
cclusttemp = permute(cclusttemp,[1 3 2]);

cellcorr.corrpossum(1,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrpossum),1);
cellcorr.corrposmean(1,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrposmean),1);
cellcorr.corrplacepossum(1,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrplacepossum),1);
cellcorr.corrplaceposmean(1,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrplaceposmean),1);
cellcorr.corrstimpossum(1,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrstimpossum),1);
cellcorr.corrstimposmean(1,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrstimposmean),1);

%% sum of correltaed place to stim cells
cellID = isplace(:,1)>0|isboth(:,1);
cclusttemp = cclust(cellID,:,:);
cclusttemp = permute(cclusttemp,[1 3 2]);

cellcorr.corrpossum(2,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrpossum),1);
cellcorr.corrposmean(2,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrposmean),1);
cellcorr.corrplacepossum(2,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrplacepossum),1);
cellcorr.corrplaceposmean(2,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrplaceposmean),1);
cellcorr.corrstimpossum(2,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrstimpossum),1);
cellcorr.corrstimposmean(2,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrstimposmean),1);

%% sum of correlated stim cells to stim cells
cellID = isstim(:,1)|isboth(:,1);
cclusttemp = cclust(cellID,:,:);
cclusttemp = permute(cclusttemp,[1 3 2]);

cellcorr.corrpossum(3,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrpossum),1);
cellcorr.corrposmean(3,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrposmean),1);
cellcorr.corrplacepossum(3,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrplacepossum),1);
cellcorr.corrplaceposmean(3,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrplaceposmean),1);
cellcorr.corrstimpossum(3,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrstimpossum),1);
cellcorr.corrstimposmean(3,~emptyCAIM) = nanmean(cclusttemp(:,int,cclustID.corrstimposmean),1);

%%


for i = find(~emptyCAIM)
    j = i+expID(1)-find(~emptyCAIM,1);
    cclust = CAIMcorr(j,mouse).cclust;
    
    isplace = (cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05) & ~(cclust(:,cclustID.nstim)>2 & cclust(:,cclustID.airpp)<=.05);
    isstim  = ~(cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05) & (cclust(:,cclustID.nstim)>2 & cclust(:,cclustID.airpp)<=.05);
    isboth  = (cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05) & (cclust(:,cclustID.nstim)>2 & cclust(:,cclustID.airpp)<=.05);
    isthere = ~isnan(cclust(:,1));
    
    cellID = isthere;
    cclusttemp = cclust(cellID,:);
    TodayCorr.corrpossum(1,i) = nanmean(cclusttemp(:,cclustID.corrpossum),1);
    TodayCorr.corrposmean(1,i) = nanmean(cclusttemp(:,cclustID.corrposmean),1);
    TodayCorr.corrplacepossum(1,i) = nanmean(cclusttemp(:,cclustID.corrplacepossum),1);
    TodayCorr.corrplaceposmean(1,i) = nanmean(cclusttemp(:,cclustID.corrplaceposmean),1);
    TodayCorr.corrstimpossum(1,i) = nanmean(cclusttemp(:,cclustID.corrstimpossum),1);
    TodayCorr.corrstimposmean(1,i) = nanmean(cclusttemp(:,cclustID.corrstimposmean),1);
    
    cellID = isplace|isboth;
    cclusttemp = cclust(cellID,:);
    TodayCorr.corrpossum(2,i) = nanmean(cclusttemp(:,cclustID.corrpossum),1);
    TodayCorr.corrposmean(2,i) = nanmean(cclusttemp(:,cclustID.corrposmean),1);
    TodayCorr.corrplacepossum(2,i) = nanmean(cclusttemp(:,cclustID.corrplacepossum),1);
    TodayCorr.corrplaceposmean(2,i) = nanmean(cclusttemp(:,cclustID.corrplaceposmean),1);
    TodayCorr.corrstimpossum(2,i) = nanmean(cclusttemp(:,cclustID.corrstimpossum),1);
    TodayCorr.corrstimposmean(2,i) = nanmean(cclusttemp(:,cclustID.corrstimposmean),1);
    
    cellID = isstim|isboth;
    cclusttemp = cclust(cellID,:);
    TodayCorr.corrpossum(3,i) = nanmean(cclusttemp(:,cclustID.corrpossum),1);
    TodayCorr.corrposmean(3,i) = nanmean(cclusttemp(:,cclustID.corrposmean),1);
    TodayCorr.corrplacepossum(3,i) = nanmean(cclusttemp(:,cclustID.corrplacepossum),1);
    TodayCorr.corrplaceposmean(3,i) = nanmean(cclusttemp(:,cclustID.corrplaceposmean),1);
    TodayCorr.corrstimpossum(3,i) = nanmean(cclusttemp(:,cclustID.corrstimpossum),1);
    TodayCorr.corrstimposmean(3,i) = nanmean(cclusttemp(:,cclustID.corrstimposmean),1);    
    
end

end
