% load('Z:\Martin Pofahl\BigFatCluster_Cues_4AP.mat')
% Airpuff = load('Z:\Martin Pofahl\Airpuff.mat');
% mouseID = Airpuff.Response.airpuff.mouseID;
%%
load('cclustID.mat','cclustID');
expID = cclustID.expID;
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
% trials to use
use = [3 4 5 ];% 

% Mice to use cells from
mice = [446 453 184 177 ];
% mice = [195 226 227];

cellID = cclust(:,expID,:);
cellID = max(cellID,[],3);
cellID = cellID == mice;
cellID = max(cellID,[],2);
% threshold of appearance
thresh = 1;

samecelltemp = samecell(use,cellID);
actcell = samecelltemp;
actcell(actcell>0) = 1;
actcell = sum(actcell,1);
% plot(actcell)
% histogram(actcell)

actcell(actcell<thresh) = 0;
actcell(actcell>0) = 1;


cclusttemp = cclust;
exclude = samecell' == 0;

for i = 1:size(cclusttemp,3)
    cclusttemp(exclude(:,i),:,i) = NaN;
end
cclusttemp = cclust(cellID,:,:);
cclusttemp = cclusttemp(actcell == 1,:,use);
samecelltemp = samecelltemp(:,actcell == 1);

%% Place cells over days

isplace = cclusttemp(:,plcvct,:)>0;% & cclusttemp(:,plcvctp,:)<=.5;

isplace1 = cclusttemp(:,plcvct,:)>0;% & cclusttemp(:,plcvctp,:)<=.5;
trans = zeros(size(isplace,3)-1,1);  %transition 
for i = 2:size(isplace,3)
    istrans = isplace1(:,1,i-1)+isplace(:,1,i);
    istrans = istrans == 2;
    trans(i-1) = sum(istrans)/sum(isplace(:,1,i)) ;
end    

isplace = sum(double(isplace(:,1,:)),3);
isplace = isplace>=1;%thresh-1;
% isplace = isplace(:,1,1)==1;

% newplace = cclusttemp(:,plcvct,end)>0 & cclusttemp(:,plcvctp,end)<=.05;
% oldplace = cclusttemp(:,plcvct,1:end-1)>0 ;%& cclusttemp(:,plcvctp,1:end-1)<=.05;
% oldplace = max(oldplace(:,1,:),[],3);
% newplace(oldplace==1) = 0;
% isplace = newplace;

plcfieldtemp = plcfield(cellID,:,:);
plcfieldtemp = plcfieldtemp(actcell == 1,:,use);
plcfieldtemp = plcfieldtemp(isplace,:,:);
b = cclusttemp(isplace,:,:);

for i = 1:size(plcfieldtemp,1)
    for j = 1:size(plcfieldtemp,3)
        plcfieldtemp(i,:,j) = (plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)))./max(plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)));
    end
end

plccenter = b(:,plcvang,1);
[~,aa] = sort(plccenter);
plccenter = plccenter(aa);
c = find(abs(plccenter)==min(abs(plccenter)));
aa = aa([c:end 1:c-1]);
% plcfieldtemp = plcfieldtemp(aa,:,:);

figure('position',[243 326 1569 645],'color',[1 1 1])
range = nanmean(plcfieldtemp,1);
range = [min(range(~isnan(range))) max(range(~isnan(range)))];
for i = 1:size(plcfieldtemp,3)
    subplot(11,size(plcfieldtemp,3),i:size(plcfieldtemp,3):9*size(plcfieldtemp,3) )
    imagesc(plcfieldtemp(:,:,i))
%     colormap(jet)
    if i == 1
%         title('All place coding cells')
    end
%     if i == size(plcfieldtemp,3)
%         hold on
%         plot([66.66 66.66],[0 size(plcfieldtemp,1)],'color',[1 0 1])
%         hold off
%     end
    ax = gca;
    ax.XTickLabel = {};
    caxis([0 1])
    subplot(11,size(plcfieldtemp,3),i+(9*size(plcfieldtemp,3)))
    imagesc(nanmean(plcfieldtemp(:,:,i),1))
%     if i == size(plcfieldtemp,3)
%         hold on
%         plot([66.66 66.66],[.5 1.5],'color',[1 0 1])
%         hold off
%     end
    caxis(range)
    ax = gca;
    ax.YTickLabel = {};
end


% stimplc = Airpuff.Response.airpuff.stimplc;

% subplot(11,size(plcfieldtemp,3),11+(9*size(plcfieldtemp,3)))
% imagesc(mean(stimplc(cell2mat(mouseID(:,2))==1,:)))
% ax = gca;
% ax.YTickLabel = {};
% ylabel('AP position')
% subplot(11,size(plcfieldtemp,3),12+(9*size(plcfieldtemp,3)))
% imagesc(mean(stimplc(cell2mat(mouseID(:,2))==2,:)))
% ax = gca;
% ax.YTickLabel = {};
% subplot(11,size(plcfieldtemp,3),13+(9*size(plcfieldtemp,3)))
% imagesc(mean(stimplc(cell2mat(mouseID(:,2))==3,:)))
% ax = gca;
% ax.YTickLabel = {};
% % printpdf('All place cells',1)

% figure
% imagesc(permute(nanmean(plcfieldtemp,1),[3 2 1]));
% imagesc(permute(nanmean(plcfieldtemp(:,:,end),1),[3 2 1]));

%% New place cells every day

isplace = cclusttemp(:,plcvct,:)>0 & cclusttemp(:,plcvctp,:)<=.05;
isplace = max(isplace(:,1,:),[],3);

plcfieldtemp = plcfield(cellID,:,:);
plcfieldtemp = plcfieldtemp(actcell == 1,:,use);
plcfieldtemp = plcfieldtemp(isplace,:,:);
b = cclusttemp(isplace,:,:);

for j = 1:size(plcfieldtemp,3)
    for i = 1:size(plcfieldtemp,1)
        plcfieldtemp(i,:,j) = (plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)))./max(plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)));
    end
end

for j = 2:size(plcfieldtemp,3)
%     newplace = cclusttemp(:,plcvct,j)>0 & cclusttemp(:,plcvctp,j)<=.05;
    oldplace = b(:,plcvct,1:j-1)>0 ;%& cclusttemp(:,plcvctp,1:end-1)<=.05;
    oldplace = max(oldplace(:,1,:),[],3);
%     newplace(oldplace==1) = 0;
    plcfieldtemp(oldplace==1,:,j) = 0;
end


plccenter = b(:,plcvang,end);
[~,aa] = sort(plccenter);
c = find(abs(plccenter)==min(abs(plccenter)));
aa = aa([c:end 1:c-1]);
aaa = aa>0;
plcfieldtemp = plcfieldtemp(aa,:,:);


figure('position',[243 526 1569 445])
title('Place coding cells')
range = nanmean(plcfieldtemp(:,:,2:end),1);
range = [min(range(:)) max(range(:))];
for i = 1:size(plcfieldtemp,3)
    subplot(10,size(plcfieldtemp,3),i:size(plcfieldtemp,3):9*size(plcfieldtemp,3) )
    imagesc(plcfieldtemp(:,:,i))
    if i == 1
        title('New place coding cells (all mice)')
    end
    if i == size(plcfieldtemp,3)
        hold on
        plot([66.66 66.66],[0 size(plcfieldtemp,2)],'color',[1 0 1])
        hold off
    end
    ax = gca;
    ax.XTickLabel = {};
    caxis([0 1])
    subplot(10,size(plcfieldtemp,3),i+(9*size(plcfieldtemp,3)))
    imagesc(nanmean(plcfieldtemp(:,:,i),1))
    if i == size(plcfieldtemp,3)
        hold on
        plot([66.66 66.66],[.5 1.5],'color',[1 0 1])
        hold off
    end
%     caxis(range)
    ax = gca;
    ax.YTickLabel = {};
end

% printpdf('new place cells',1)

% figure
% a = permute(nanmean(plcfieldtemp(:,:,1:end),1),[3 2 1]);
% subplot(1,2,1)
% imagesc(a(4:end,:))
% subplot(1,2,2)
% plot(a(4:end,:)')
%% runperc cells

a =  cclusttemp(:,runperc,:)>.8 & cclusttemp(:,runnerp,:)<.05 & cclusttemp(:,meanf,:)>.5;
% a = cclusttemp(:,meanf,:)>1;
b = cclusttemp(max(a(:,1,:),[],3),:,:);

c = b(:,runperc,:);
c = reshape(c,size(c,1),size(c,3));
% c = c(~isnan(c));
% histogram(c,100)

cerr = nanstd(c);
cerr = cerr/sqrt(size(c,2));
errorbar(1:size(c,2),nanmean(c),cerr)


