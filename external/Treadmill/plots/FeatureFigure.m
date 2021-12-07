load('Z:\Martin Pofahl\BigFatCluster_Cues_4AP.mat')
load('cclustID.mat')

%%
% trials to use
use = [6:11];% 
% Mice to use cells from
mice = [195 226 227 ];
% mice = [103 194 195 226 227];
% mice = [194 195 224 226 227 234];
% mice = [103 226 227];
cellID = cclust(:,cclustID.expID,:);
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

%% open figure
figure('color',[1 1 1],...
    'position',[300 0 2*[420 594]],...
    'visible','on',...
    'PaperUnits','centimeters',...
    'Units','centimeters')

%% number of feature cells per day


isplaceday = (cclusttemp(:,cclustID.plcvct,:)>0 & cclusttemp(:,cclustID.plcvctp,:)<=.05);
isplaceday = permute(isplaceday,[1 3 2]);
isplaceday(cclusttemp(:,cclustID.followcell,1)==0,:) = 0;
isplaceday = [isplaceday(:,1)|isplaceday(:,2)|isplaceday(:,3) isplaceday(:,4)|isplaceday(:,5)|isplaceday(:,6)];

isstimday = (cclusttemp(:,cclustID.nstim,:)>4 & cclusttemp(:,cclustID.airpp,:)<=.05);
isstimday = permute(isstimday,[1 3 2]);
isstimday(cclusttemp(:,cclustID.followcell,1)==0,:) = 0;
isstimday = [isstimday(:,1)|isstimday(:,2)|isstimday(:,3) isstimday(:,4)|isstimday(:,5)|isstimday(:,6)];

isbothday = isplaceday & isstimday;
isplaceday(isbothday==1) = 0;
isstimday(isbothday==1) = 0;



% subplot(5,1,1)
% p = bar([sum(isstimday,1)' sum(isbothday,1)' sum(isplaceday,1)']);
% piecol = {[1 0 1]     ,[1 1 0]           ,[0 1 1]    ,[.5 .5 .5]};
% for i = 1:length(p)
%     pp = p(i);
%     pp.FaceColor  = piecol{i};
% end
% title('Number of feature cells (day code)')
% ax = gca;
% ax.XTickLabel = {'Cue1' 'Cue2' 'Cue3' 'AP1' 'AP2' 'AP3' 'AP fixed'};
% xtickangle(45)


%% feature cells over days 

isplace = (cclusttemp(:,cclustID.plcvct,:)>0 & cclusttemp(:,cclustID.plcvctp,:)<=.05);
isplace = sum(isplace,3);
isplace(isplace>1) = 1;
noplace = ~isplace;
isplace = isplace & cclusttemp(:,cclustID.followcell,1) == 1;
noplace = noplace & cclusttemp(:,cclustID.followcell,1) == 1;

isstim = (cclusttemp(:,cclustID.nstim,4)>4 & cclusttemp(:,cclustID.airpp,4)<=.05)|(cclusttemp(:,cclustID.nstim,5)>4 & cclusttemp(:,cclustID.airpp,5)<=.05)|(cclusttemp(:,cclustID.nstim,6)>4 & cclusttemp(:,cclustID.airpp,6)<=.05);
nostim = ~isstim;
isstim = isstim & cclusttemp(:,cclustID.followcell,1) == 1;
nostim = nostim & cclusttemp(:,cclustID.followcell,1) == 1;

isboth = isplace&isstim;
isnone = nostim&noplace;
isplace(isboth) = 0;
isstim(isboth) = 0;


%% feature change over days


b = [find(isstim); find(isboth); find(isplace)];
subplot(5,1,1)
% bar([sum(isstim) sum(isboth) sum(isplace)],'stacked')
hold on
for i = 1:sum(isstim)+sum(isboth)+sum(isplace)
%     plot([1 7],[i i],'color',[.5 .5 .5])
    c = find(isstimday(b(i),:));
    scatter(c,i*ones(length(c),1),10,'filled','MarkerFaceColor',[1 0 1])   
    c = find(isbothday(b(i),:));
    scatter(c,i*ones(length(c),1),10,'filled','MarkerFaceColor',[1 1 0])
    c = find(isplaceday(b(i),:));
    scatter(c,i*ones(length(c),1),10,'filled','MarkerFaceColor',[0 1 1])
end
axis tight
xlim([0.5 3.5])
% xlim([0.5 7.5])
ax = gca;
ax.XTick = 1:3;
ax.XTickLabel = {'Cues' 'AP' 'AP fixed'};
% ax.XTickLabel = {'Cue1' 'Cue2' 'Cue3' 'AP1' 'AP2' 'AP3' 'AP fixed'};
xtickangle(45)
title('Cell feature appearance')
%%

subplot(5,1,2)

feature = cclustID.netprob; % cclustID.meanf cclustID.netprob cclustID.netperc

a = cclusttemp(isstim,feature,:);%/sum(isstim);
a(a==0) = NaN;
a = permute(a,[1 3 2]);
a = [nanmean(a(:,1:3),2) nanmean(a(:,4:6),2) a(:,7)];
aerr = nanstd(a,1)/sqrt(sum(isstim));
a = nanmean(a,1);


b = cclusttemp(isboth,feature,:);
b(b==0) = NaN;
b = permute(b,[1 3 2]);
b = [nanmean(b(:,1:3),2) nanmean(b(:,4:6),2) b(:,7)];
berr = nanstd(b,1)/sqrt(sum(isboth));
b = nanmean(b,1);%/sum(isplace);


c = cclusttemp(isplace,feature,:);
c(c==0) = NaN;
c = permute(c,[1 3 2]);
c = [nanmean(c(:,1:3),2) nanmean(c(:,4:6),2) c(:,7)];
cerr = nanstd(c,1)/sqrt(sum(isboth));
c = nanmean(c,1);%/sum(isplace);

d = cclusttemp(isnone,feature,:);
d(d==0) = NaN;
d = permute(d,[1 3 2]);
d = [nanmean(d(:,1:3),2) nanmean(d(:,4:6),2) d(:,7)];
derr = nanstd(d,1)/sqrt(sum(isboth));
d = nanmean(d,1);%/sum(isplace);

p = bar([a' b' c' d']);
hold on
errorbar((1:length(a))-.27,a,aerr,'.','Color',[0 0 0])
errorbar((1:length(a))-.1,b,berr,'.','Color',[0 0 0])
errorbar((1:length(a))+.1,c,cerr,'.','Color',[0 0 0])
errorbar((1:length(a))+.27,d,derr,'.','Color',[0 0 0])
piecol = {[1 0 1]     ,[1 1 0]           ,[0 1 1]    ,[.5 .5 .5]};
for i = 1:length(p)
    pp = p(i);
    pp.FaceColor  = piecol{i};
end

title('Network Probability')
ax = gca;
ax.XTickLabel = {'Cues' 'AP' 'AP fixed'};
ylabel('fraction of participated net-ev')
% ax.XTickLabel = {'Cue1' 'Cue2' 'Cue3' 'AP1' 'AP2' 'AP3' 'AP fixed'};
xtickangle(45)
legend(p,{'AP' 'Both' 'Place' 'None'},'Location','best','Orientation','horizontal')
legend('boxoff')
% histogram(cclusttemp(isstim,cclustID.netprob,:))
% hold on
% histogram(cclusttemp(isplace,cclustID.netprob,:))
% hold off
%%
feature = cclustID.meanf; % cclustID.meanf cclustID.netprob cclustID.netperc
subplot(5,1,3)

a = cclusttemp(isstim,feature,:);%/sum(isstim);
a(a==0) = NaN;
a = permute(a,[1 3 2]);
a = [nanmean(a(:,1:3),2) nanmean(a(:,4:6),2) a(:,7)];
aerr = nanstd(a,1)/sqrt(sum(isstim));
a = nanmean(a,1);


b = cclusttemp(isboth,feature,:);
b(b==0) = NaN;
b = permute(b,[1 3 2]);
b = [nanmean(b(:,1:3),2) nanmean(b(:,4:6),2) b(:,7)];
berr = nanstd(b,1)/sqrt(sum(isboth));
b = nanmean(b,1);%/sum(isplace);


c = cclusttemp(isplace,feature,:);
c(c==0) = NaN;
c = permute(c,[1 3 2]);
c = [nanmean(c(:,1:3),2) nanmean(c(:,4:6),2) c(:,7)];
cerr = nanstd(c,1)/sqrt(sum(isboth));
c = nanmean(c,1);%/sum(isplace);

d = cclusttemp(isnone,feature,:);
d(d==0) = NaN;
d = permute(d,[1 3 2]);
d = [nanmean(d(:,1:3),2) nanmean(d(:,4:6),2) d(:,7)];
derr = nanstd(d,1)/sqrt(sum(isboth));
d = nanmean(d,1);%/sum(isplace);

p = bar([a' b' c' d']);
hold on
errorbar((1:length(a))-.27,a,aerr,'.','Color',[0 0 0])
errorbar((1:length(a))-.1,b,berr,'.','Color',[0 0 0])
errorbar((1:length(a))+.1,c,cerr,'.','Color',[0 0 0])
errorbar((1:length(a))+.27,d,derr,'.','Color',[0 0 0])
piecol = {[1 0 1]     ,[1 1 0]           ,[0 1 1]    ,[.5 .5 .5]};
for i = 1:length(p)
    pp = p(i);
    pp.FaceColor  = piecol{i};
end
title('Mean frequency')
ax = gca;
% ax.XTickLabel = [];
ax.YLim = [0 20];
ax.XTickLabel = {'Cues' 'AP' 'AP fixed'};
ylabel('Events min^-^1')
% ax.XTickLabel = {'Cue1' 'Cue2' 'Cue3' 'AP1' 'AP2' 'AP3' 'AP fixed'};
xtickangle(45)
%%
feature = cclustID.plcvct; % cclustID.meanf cclustID.netprob cclustID.netperc
subplot(5,1,4)
a = cclusttemp(isstim,feature,:);%/sum(isstim);
a(a==0) = NaN;
a = permute(a,[1 3 2]);
a = [nanmean(a(:,1:3),2) nanmean(a(:,4:6),2) a(:,7)];
aerr = nanstd(a,1)/sqrt(sum(isstim));
a = nanmean(a,1);


b = cclusttemp(isboth,feature,:);
b(b==0) = NaN;
b = permute(b,[1 3 2]);
b = [nanmean(b(:,1:3),2) nanmean(b(:,4:6),2) b(:,7)];
berr = nanstd(b,1)/sqrt(sum(isboth));
b = nanmean(b,1);%/sum(isplace);


c = cclusttemp(isplace,feature,:);
c(c==0) = NaN;
c = permute(c,[1 3 2]);
c = [nanmean(c(:,1:3),2) nanmean(c(:,4:6),2) c(:,7)];
cerr = nanstd(c,1)/sqrt(sum(isboth));
c = nanmean(c,1);%/sum(isplace);

d = cclusttemp(isnone,feature,:);
d(d==0) = NaN;
d = permute(d,[1 3 2]);
d = [nanmean(d(:,1:3),2) nanmean(d(:,4:6),2) d(:,7)];
derr = nanstd(d,1)/sqrt(sum(isboth));
d = nanmean(d,1);%/sum(isplace);

p = bar([a' b' c' d']);

hold on
errorbar((1:length(a))-.27,a,aerr,'.','Color',[0 0 0])
errorbar((1:length(a))-.1,b,berr,'.','Color',[0 0 0])
errorbar((1:length(a))+.1,c,cerr,'.','Color',[0 0 0])
errorbar((1:length(a))+.27,d,derr,'.','Color',[0 0 0])
piecol = {[1 0 1]     ,[1 1 0]           ,[0 1 1]    ,[.5 .5 .5]};
for i = 1:length(p)
    pp = p(i);
    pp.FaceColor  = piecol{i};
end
title('Place vector length')
ax = gca;
% ax.XTickLabel = [];
ax.XTickLabel = {'Cues' 'AP' 'AP fixed'};
ylabel('Normalized vector length')
% ax.XTickLabel = {'Cue1' 'Cue2' 'Cue3' 'AP1' 'AP2' 'AP3' 'AP fixed'};
xtickangle(45)
%% transition of cells

% isplace = (cclusttemp(:,cclustID.plcvct,:)>0 & cclusttemp(:,cclustID.plcvctp,:)<=.05);
% isplace = permute(isplace,[1 3 2]);
% isstim = (cclusttemp(:,cclustID.nstim,:)>4 & cclusttemp(:,cclustID.airpp,:)<=.05);
% isstim = permute(isstim,[1 3 2]);
% isboth = isplace&isstim;

no2bo = zeros(size(isplaceday,1),1);
pl2bo = zeros(size(isplaceday,1),1);
ap2bo = zeros(size(isplaceday,1),1);
bo2no = zeros(size(isplaceday,1),1);
bo2pl = zeros(size(isplaceday,1),1);
bo2ap = zeros(size(isplaceday,1),1);
ap2ap = zeros(size(isplaceday,1),1);
pl2pl = zeros(size(isplaceday,1),1);
no2no = zeros(size(isplaceday,1),1);
pl2ap = zeros(size(isplaceday,1),1);
ap2pl = zeros(size(isplaceday,1),1);
no2ap = zeros(size(isplaceday,1),1);
ap2no = zeros(size(isplaceday,1),1);
no2pl = zeros(size(isplaceday,1),1);
pl2no = zeros(size(isplaceday,1),1);

for i = 1:size(isplaceday,1)
    for j = 1:size(isplaceday,2)-1
        if (isplaceday(i,j) == 0 && isstimday(i,j) == 0) && isbothday(i,j+1) == 1
            no2bo(i) = no2bo(i) + 1;
        elseif isplaceday(i,j) == 1 && isbothday(i,j+1) == 1
            pl2bo(i) = pl2bo(i) + 1;
        elseif isstimday(i,j) == 1 && isbothday(i,j+1) == 1
            ap2bo(i) = ap2bo(i) + 1;
        elseif isbothday(i,j) == 1 && (isstimday(i,j+1) == 0 && isplaceday(i,j+1) == 0)
            bo2no(i) = bo2no(i) + 1;
        elseif isbothday(i,j) == 1 && isplaceday(i,j+1) == 1
            bo2pl(i) = bo2pl(i) + 1;
        elseif isbothday(i,j) == 1 && isstimday(i,j+1) == 1
            bo2ap(i) = bo2ap(i) + 1;
        elseif isstimday(i,j) == 1 && isstimday(i,j+1) ==1
            ap2ap(i) = ap2ap(i) + 1;
        elseif isplaceday(i,j) == 1 && isplaceday(i,j+1)== 1
            pl2pl(i) = pl2pl(i) +1;
        elseif (isplaceday(i,j) == 0 && isstimday(i,j) == 0) && (isplaceday(i,j+1) == 0 && isstimday(i,j+1) == 0)
            no2no(i) = no2no(i) +1;
        elseif isplaceday(i,j)== 1 && isstimday(i,j+1)== 1
            pl2ap(i) = pl2ap(i) + 1;
        elseif isstimday(i,j)== 1 && isplaceday(i,j+1) == 1
            ap2pl(i) = ap2pl(i) + 1;
        elseif (isplaceday(i,j) == 0 && isstimday(i,j) == 0) && isstimday(i,j+1) == 1 
            no2ap(i) = no2ap(i) + 1;
        elseif isstimday(i,j) == 1 && (isplaceday(i,j+1) == 0 && isstimday(i,j+1) == 0)
            ap2no(i) = ap2no(i) + 1;
        elseif (isstimday(i,j) == 0 && isplaceday(i,j) == 0) && isplaceday(i,j+1) == 1
            no2pl(i) = no2pl(i) + 1;
        elseif isplaceday(i,j) == 1 && (isplaceday(i,j+1) == 0 && isstimday(i,j+1) == 0)
            pl2no(i) = pl2no(i) + 1;
        end
    end
end

a = [sum(no2bo) sum(pl2bo) sum(ap2bo) sum(bo2no) sum(bo2pl) sum(bo2ap) sum(ap2ap) sum(pl2pl) sum(pl2ap) sum(ap2pl) sum(no2ap) sum(ap2no) sum(no2pl) sum(pl2no)];
subplot(5,1,5)
bar(a)
ax = gca;
ax.XTickLabel = {'no2bo' 'pl2bo' 'ap2bo' 'bo2no' 'bo2pl' 'bo2ap' 'ap2ap' 'pl2pl' 'pl2ap' 'ap2pl' 'no2ap' 'ap2no' 'no2pl' 'pl2no'};
xtickangle(45)
title('Change of feature')
grid on


%%
printpdf('Cell features',1)
