%% CellIDplot
Cues = load('Z:\Martin Pofahl\Cues.mat');
Baseline = load('Z:\Martin Pofahl\Baseline.mat');
Airpuff = load('Z:\Martin Pofahl\Airpuff.mat');
load('Z:\Martin Pofahl\BigFatCluster.mat');
load('Z:\Martin Pofahl\BigFatCluster_Bs_Cue.mat');
load('cclustID.mat')
%%
experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};
mouse = {'M103' 'M155' 'M158' 'M194' 'M195' 'M224' 'M226' 'M227' 'M229' 'M234'}; 
plotcol = [0 .3 0;
    0 .6 0;
    0 .9 0;
    .3 0 0;
    .6 0 0;
    .9 0 0;
    0 0 .3;
    0 0 .6;
    0 0 .9];

trial = 3;
ftsz = 8;
lnwd = 1;
close all

%%
figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]],...
        'renderer','painters',...
        'visible','on')
%% mean firing frequency of active cells

figure

x = [1 2 4 5 7 8 10 11 13 14 16 17];
mouseID = Baseline.Fire.mouseID;
yin = Baseline.Fire.meanfire(:,2:3);
y = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
yerr = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,2)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,2))==i,:),1);
    yerr(i,:) = nanstd(yin(cell2mat(mouseID(:,2))==i,:),1);
end

y = y(3:5,:);
y = reshape(y',size(y,1)*size(y,2),1);
yerr = yerr(3:5,:);
yerr = reshape(yerr',size(yerr,1)*size(yerr,2),1);

mouseID = Cues.Fire.mouseID;
yin =  Cues.Fire.meanfire(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
y1err = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,2)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,2))==i,:),1);
    y1err(i,:) = nanstd(yin(cell2mat(mouseID(:,2))==i,:),1);
end

y1 = reshape(y1',size(y1,1)*size(y1,2),1);
y1err = reshape(y1err',size(y1err,1)*size(y1err,2),1);

y = [y; y1];
yerr = [yerr; y1err];
yerr = yerr./sqrt(max(cell2mat(mouseID(:,4))));

% reps = size(y,1);
% [~,~,stats] = anova2(ystat,reps,'off');
% c = multcompare(stats,'Display','off');
% c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

b = bar(x,y);
hold on

errorbar(x,y,yerr,'.',...
            'Color',[0 0 0],...
            'LineWidth',lnwd)

b.FaceColor = 'flat';
b.LineStyle = 'none';
for i = 1:length(y)/2
    b.CData(2*i-1,:) = [.6 .6 .6];
    b.CData(2*i,:) = [.2 .2 .2];
end

ax = gca;
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XTick = 1.5 :3: 16.5;
ax.XTickLabel = {'Base1' 'Base2' 'Base3' 'Cues1' 'Cues2' 'Cues3'};

title('Mean firing frequency')


% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\figure3a'); 

%% stacked bar graph with cell identities 

figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]],...
        'renderer','painters',...
        'visible','on')
    
% subplot(3,2,1)
use = [3 4 5 6 7 8];% 
mice = [103 155 158 194 195 224 226 227 234];% Mice to use cells from
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

actcell(actcell<thresh) = 0;
actcell(actcell>0) = 1;


% exclude = samecell' == 0;
% for i = 1:size(cclusttemp,3)
%     cclusttemp(exclude(:,i),:,i) = NaN;
% end
% cclusttemp = cclust(cellID,:,:);
% cclusttemp = cclusttemp(actcell == 1,:,:);
% samecelltemp = samecelltemp(:,actcell == 1);

% %%
a = zeros(length(use),8);
for i = 1:length(use)

    cclusttemp = cclust(:,:,use(i));
    cclusttemp = cclusttemp(cclusttemp(:,cclustID.expID)>0,:);
    isnet = cclusttemp(:,cclustID.netprob)>0;
    isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
    isstim = cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05;
    
    % network, no place, stim
%     a(i,1) = sum(isnet &  ~isplace & isstim)/length(cclusttemp);
    % no network, no place, stim
%     a(i,2) = sum(~isnet &  ~isplace & isstim)/length(cclusttemp);
    % network, place, stim
%     a(i,3) = sum(isnet &  isplace & isstim)/length(cclusttemp);
    % no network, place, stim
%     a(i,4) = sum(~isnet &  isplace & isstim)/length(cclusttemp);
    % network, place, no stim
    a(i,5) = sum(isnet &  isplace & ~isstim)/length(cclusttemp);
    % no network, place, no stim
    a(i,6) = sum(~isnet &  isplace & ~isstim)/length(cclusttemp);
    % network, no place, no stim
%     a(i,7) = sum(isnet  &  ~isplace & ~isstim)/length(cclusttemp);
%     no network, no place, no stim
%     a(i,8) = sum(~isnet &  ~isplace & ~isstim)/length(cclusttemp);
end
a = (100*a);



p = bar(a,'stacked');
labels = {'net & stim','stim','net & plc & stim','plc & stim','net & plc','plc','net','non'};
piecol = {[1 0 1]     ,[.5 0 .5],[1 1 0]        ,[.5 .5 0]   ,[0 1 1]    ,[0 .5 .5],[0 1 0],[.5 .5 .5]};

for i = 1:length(a)
    pp = p(i);
    pp.FaceColor  = piecol{i};
end
box('off')

ax = gca;
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
% ax.XTick = [1.5 4.5];
ax.XTickLabel = {'Base1' 'Base2' 'Base3' 'Cues1' 'Cues2' 'Cues3'};
ylabel('Fraction of cells %','FontSize',ftsz)
% ylim([0 .2])
% xlim([1.5 2.5])

title('PCs')
%% bar graph with place vector length 
subplot(3,2,2)

x = [1 2 4 5 7 8 10 11 13 14 16 17];

CAIM = CAIM(:,:);
plcvct = zeros(8,size(CAIM,2));
% plcerr = zeros(8,size(CAIM,2));
for i = 3:8%1:size(CAIM,1)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];
    netpos = [];
    for j = 1:length(fullCAIM)
        k = fullCAIM(j);      
        cclusttemp = CAIM(i,k).cclust;
        isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=0.05;
        isnet = cclusttemp(:,cclustID.netprob)>0;
        plcvct(i,j,1) = nanmean(cclusttemp(isplace & isnet,cclustID.plcvct));
        plcvct(i,j,2) = nanmean(cclusttemp(isplace & ~isnet,cclustID.plcvct));
%         plcerr(i,j,1) = nanstd(cclusttemp(isplace & isnet,cclustID.plcvct));
%         plcerr(i,j,2) = nanstd(cclusttemp(isplace & ~isnet,cclustID.plcvct));
    end
end

% plcerr(plcerr == 0) = NaN;
% plcerr = plcerr./sqrt(size(CAIM,2));
plcvct(plcvct == 0) = NaN; 
y = [nanmean(plcvct(3:end,:,1),2) nanmean(plcvct(3:end,:,2),2)];
yerr = [nanstd(plcvct(3:end,:,1),[],2) nanstd(plcvct(3:end,:,2),[],2)];
y = reshape(y',1,size(y,1)*2);
yerr = reshape(yerr',1,size(yerr,1)*2);

b = bar(x,y);
hold on

errorbar(x,y,yerr,'.',...
            'Color',[0 0 0],...
            'LineWidth',lnwd)

b.FaceColor = 'flat';
b.LineStyle = 'none';
for i = 1:length(y)/2
    b.CData(2*i-1,:) = [0 1 1];
    b.CData(2*i,:) = [0 .5 .5];
end

ax = gca;
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XTick = 1.5 :3: 16.5;
ax.XTickLabel = {'Base1' 'Base2' 'Base3' 'Cues1' 'Cues2' 'Cues3'};

hold off
title('Place Vector length')

%% single cells speed correlation
figure('color',[1 1 1],...
        'position',[276   363   944   631],...[500 50 1.5*[420 594]],...
        'renderer','painters',...
        'visible','on')
 aa = [];
for i = 3:8%1:size(CAIM,1)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];

    for j = 1:length(fullCAIM)
        
        k = fullCAIM(j); 

        speedcorr = CAIM(i,k).speedcorr;
        if ~isempty(speedcorr)
            cellID = speedcorr.cellID;
            trace = speedcorr.trace;
            isspeed = cellID(:,1)>.9 & cellID(:,2)<.05;
            isneg = cellID(:,1)<-.9 & cellID(:,2)<.05;

    %         subplot(6,18,(i-3)*18+(2*k-1))  
            subplot(6,9,(i-3)*9+k)  

            plot(trace(isspeed,:)','color',[.5 .5 .5])
            hold on
            plot(mean(trace(isspeed,:)),'color',[0 0 1],'linewidth',2)
    %         ax1 = gca;
    %         subplot(6,18,(i-3)*18+(2*k))  
%             plot(trace(isneg,:)','color',[.5 .5 .5])
%             hold on
%             plot(mean(trace(isneg,:)),'color',[1 0 0],'linewidth',2)
    %         hold off
    %         ax2 = gca;
    %         ylim = [min([ax1.YLim(1) ax2.YLim(1)]) max([ax1.YLim(2) ax2.YLim(2)])];
    %         ax2.YLim = ylim;
    %         ax1.YLim = ylim;
            
            cclusttemp = CAIM(i,k).cclust;
            isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
            aa = [aa;[sum(isplace & isspeed) sum(isplace & isneg)]];
        end
    end    
end

%% basic properties of amplitude stats
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

num = 3:8;
SumStat = zeros(length(num),size(CAIM,2),6);
a = zeros(2,1);
b = zeros(10,2);

for i = num
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];    
    C = [];
    % 1:2 variance 3:4 skewness 5:6 kurtosis
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        c = CAIM(i,k).fireprop.BsStat(:,1:2);
        C = [C; c];
        SumStat(i-num(1)+1,k,:) = CAIM(i,k).fireprop.SumStat;
        a(:,1) = SumStat(i-num(1)+1,k,1:2);
        subplot(6,11,(i-num(1))*11+k) 
        bar(a)
%         boxplot(c);
        ax = gca;
        ax.XTickLabel = {'rest' 'run'};
        if i == 3;title(mouse(k));end
%         
    end
    
    subplot(6,11,(i-num(1))*11+11) 
    b(:,:) = SumStat(i-num(1)+1,:,3:4);
    boxplot(b);
%     boxplot(C);
    ax = gca;
    ax.XTickLabel = {'rest' 'run'};
    hold off
    if i == 3;title('pool');end
end


% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\Stats\Skewness ind');
%% individual cell correlation
num = 4:5;
cclust = [];
for i = num
    %%
    fullCAIM = 1:size(CAIMcorr,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIMcorr(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    
    %%
    for j = 1:length(fullCAIM)
        k = fullCAIM(j);
        temp = CAIM(i,k).cclust;
        temp = temp(~isnan(temp(:,1)),:);
        cclust = [cclust; temp];       
        
    end
    if ~isempty(cclust)
        cellcorrplot(cclust)
%         printpdf([experiment{i} ' correlations'])
    end
end
    