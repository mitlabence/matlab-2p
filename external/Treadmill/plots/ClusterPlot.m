%% Cluster Plots 

Cues = load('Z:\Martin Pofahl\Cues.mat');
Baseline = load('Z:\Martin Pofahl\Baseline.mat');
Airpuff = load('Z:\Martin Pofahl\Airpuff.mat');
load('Z:\Martin Pofahl\BigFatCluster.mat');
% load('Z:\Martin Pofahl\BigFatCluster_Bs_Cue.mat');
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




%% Cumsum cluster size

subplot(2,1,1)
hold on

for i = 3:8%1:size(CAIM,1)
    %%
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    clustsize = [];

    for j = 1:length(fullCAIM)
        k = fullCAIM(j);      
        if ~isempty(CAIM(i,k).network.clustsize)
            clustsize = [clustsize; histcounts(CAIM(i,k).network.clustsize,1:41)];
        end
    end
    
    plot(1:40,cumsum(sum(clustsize)/sum(sum(clustsize))),...
    'LineWidth',lnwd,...
    'color',plotcol(i-2,:));
end
legend(experiment(3:8),'Location','best')
title('Cluster size')
%% Cumsum PC cluster size

subplot(2,1,2)
hold on

for i = 3:8%1:size(CAIM,1)
    %%
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    clustsizePC = [];

    for j = 1:length(fullCAIM)
        k = fullCAIM(j);    
        if ~isempty(CAIM(i,k).network.clustsize)
            clustsizePC = [clustsizePC; CAIM(i,k).network.clustsize(2,:)];
        end
    end
    
    plot(1:40,cumsum(sum(clustsizePC)/sum(sum(clustsizePC))),...
    'LineWidth',lnwd,...
    'color',plotcol(i-2,:));
end
legend(experiment(3:8),'Location','best')
title('PC Cluster size')

%%
% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\figure3b - no running networks'); 


%% Featuress

clustsize = [];
clustsizePC = [];   
clustsizenone = [];
clustsizeSpos = [];
clustsizeSposPC = [];

clustfreq = [];
% cluststim = [];
clustplace = [];
clustpos = [];
clustConcept = [];
clustsig = [];
ClustCell = [];

netnorm = [];
netnormair = [];
meanr = [];
netnormplace = [];
num = 3:5;

for i = num
    %%
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    %%

%     clustsize = [];
%     clustsizePC = [];   
%     clustsizenone = [];
%     clustsizeSpos = [];
%     clustsizeSposPC = [];
%     
%     clustfreq = [];
% %     cluststim = [];
%     clustplace = [];
%     clustpos = [];
%     clustConcept = [];
%     clustsig = [];
%     ClustCell = [];
%     
%     netnorm = [];
%     netnormair = [];
%     netnormplace = [];
    %%
    for j = 1:length(fullCAIM)
        k = fullCAIM(j);
             
        if ~isempty(CAIM(i,k).network.clustsize)
            %%
            network =  CAIM(i,k).network;
            clustID = network.clustID; 
            isCsig      = clustID.clustsig(2,:) == 1;
            isCpos      = clustID.ispos  >  0 & clustID.isplace == 0 & isCsig;
            isCposplace = clustID.ispos  >  0 & clustID.isplace  > 0 & isCsig;
            isCplace    = clustID.ispos  == 0 & clustID.isplace  > 0 & isCsig;
            isCnone     = clustID.ispos  == 0 & clustID.isplace == 0 & isCsig;
                      
            clustsize = [clustsize network.clustsize(isCsig)];
            clustsizePC = [clustsizePC network.clustsize(isCplace)];
%             clustsizeAP = [clustsizeAP; network.clustsize(isstim)];
            clustsizenone = [clustsizenone network.clustsize(isCnone)];
            clustsizeSpos = [clustsizeSpos network.clustsize(isCpos)];
%             clustsizeSneg = [clustsizeSneg; network.clustsize(isneg)];
            clustsizeSposPC = [clustsizeSposPC network.clustsize(isCposplace)];
%             clustsizeSnegPC = [clustsizeSnegPC; network.clustsize(isnegplace)];

            clustfreq = [clustfreq network.clustfreq(isCsig)];
%             cluststim = [cluststim network.cluststim];
            clustplace = [clustplace network.clustplace(isCsig)];
            clustpos = [clustpos network.clustpos(isCsig)];
            clustConcept = [clustConcept; [sum(isCplace) sum(isCposplace) sum(isCpos) sum(isCnone)]];
            clustsig = [clustsig; [sum(isCsig) sum(isCsig == 0) sum(isCsig & isCplace) sum(isCsig==0 & isCplace) sum(isCsig & isCposplace) sum(isCsig==0 & isCposplace) sum(isCsig & isCpos) sum(isCsig==0 & isCpos)]];
            
            %%
            cclusttemp = CAIM(i,k).cclust;
            isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
            isspeed = cclusttemp(:,cclustID.speedcorr)>.8 & cclusttemp(:,cclustID.speedcorrp1)<=.05;
            b = size(cclusttemp,1);
            a = histcounts(cclusttemp(:,cclustID.clust),1:length(isCsig)+1);
            b(2) = sum(a(isCsig));
            a = histcounts(cclusttemp(isplace,cclustID.clust),1:length(isCsig)+1);
            b(3) = sum(isplace); b(4) = sum(a(isCsig));
            a = histcounts(cclusttemp(isspeed,cclustID.clust),1:length(isCsig)+1);
            b(5) = sum(isspeed); b(6) = sum(a(isCsig));
            a = histcounts(cclusttemp(~isplace & ~isspeed,cclustID.clust),1:length(isCsig)+1);
            b(7) = sum(~isplace & ~isspeed); b(8) = sum(a(isCsig));
            ClustCell = [ClustCell;b];
            
            %%
            netcorr = CAIM(i,k).network.netcorr;
            a = netcorr;
            a = tril(a,-1);
            b = nanmean(a(:));
            a = netcorr(isplace,isplace);
            a = tril(a,-1);
            b(2) = nanmean(a(:));
            a = netcorr(isspeed,isspeed);
            a = tril(a,-1);
            b(3) = nanmean(a(:));
            a = netcorr(~isspeed & ~isplace,~isspeed & ~isplace);
            a = tril(a,-1);
            b(4) = nanmean(a(:));
            meanr = [meanr; b]; 
        end
    end
    
    %%
    figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]],...
        'renderer','painters',...
        'visible','on')

    subplot(3,2,1)
    hold on
    plot(1:30,cumsum(histcounts(clustsize,1:31))/sum(histcounts(clustsize,1:31)),'color',[0 0 0])
    plot(1:30,cumsum(histcounts(clustsizePC,1:31))/sum(histcounts(clustsizePC,1:31)),'color',[0 1 1])
    plot(1:30,cumsum(histcounts(clustsizeSpos,1:31))/sum(histcounts(clustsizeSpos,1:31)),'color',[1 .5 0])
    plot(1:30,cumsum(histcounts(clustsizeSposPC,1:31))/sum(histcounts(clustsizeSposPC,1:31)),'color',[1 1 0])
    plot(1:30,cumsum(histcounts(clustsizenone,1:31))/sum(histcounts(clustsizenone,1:31)),'color',[.5 .5 .5])
    title('Cluster size')
    ylabel('Clusters count')
    xlabel('Cells per cluster')
    hold off
    subplot(3,2,2)
    bar(.1:.1:1,histcounts(clustfreq,0.1:.1:1.1))
    title('Cluster freq')
    ylabel('Clusters count')
    xlabel('% of net events')   
    subplot(3,2,3)
    bar(1:10,histcounts(clustplace,1:11))
    title('Place cells in cluster')
    ylabel('Clusters count')
    xlabel('PCs per cluster')
    subplot(3,2,4)
    bar(1:10,histcounts(clustpos,1:11))
    title('Speed cells in Cluster')
    ylabel('Clusters count')
    xlabel('Speed Cells per cluster')
    
    subplot(3,2,5)
    labels = { 'plc' , 'sp & plc', 'sp' ,'none'};   
    piecol = {[0 1 1],[1 1 0]       ,[1 .5 0] ,[.5 .5 .5]};
%             isstim   isboth  isCplace isCnone    isCpos        isneg   isCposplace isnegplace ispureplace

%     labels = {'stim' ,  'both',  'plc','none',   'sppos',    'spneg','sppos & plc','spneg & plc','pure plc'};   
%     piecol = {[1 0 1],[1 1 0],[0 1 1],[.5 .5 .5],[.6 .35 0],[.35 .6 0],[1 0 0],[0 1 0],[0 0 1]};
% %             isstim   isboth  isCplace isCnone    isCpos        isneg   isCposplace isnegplace ispureplace    
    
%     figure
    clustConcept1 = sum(clustConcept,1);  
    piecol = piecol(clustConcept1>0);
    labels = labels(clustConcept1>0);
    a = clustConcept1(clustConcept1>0);
    p = pie(a,labels);
    for jj = 1:length(a)
        pp = p(2*jj-1);
        pp.FaceColor  = piecol{jj};
    end
    
    subplot(3,2,6)
    a = sum(clustsig(:,1:2:end));
    a(2,:) = sum(clustsig(:,2:2:end));
    a = a./(sum(a));
    bar(a','stacked')
    ax = gca;
    ax.XTickLabel = { 'total','plc' , 'both', 'sp' }; 
    
end

a = cumsum(histcounts(clustsize,1:31))/sum(histcounts(clustsize,1:31));
clustmed(1) = find(a>.5,1);
a = cumsum(histcounts(clustsizePC,1:31))/sum(histcounts(clustsizePC,1:31));
clustmed(2) = find(a>.5,1);
a = cumsum(histcounts(clustsizeSpos,1:31))/sum(histcounts(clustsizeSpos,1:31));
clustmed(3) = find(a>.5,1);
a = cumsum(histcounts(clustsizeSposPC,1:31))/sum(histcounts(clustsizeSposPC,1:31));
clustmed(4) = find(a>.5,1);
a = cumsum(histcounts(clustsizenone,1:31))/sum(histcounts(clustsizenone,1:31));
clustmed(5) = find(a>.5,1);

meanr = nanmean(meanr,1);

ClustCell = sum(ClustCell);
ClustCell = ClustCell(2:2:end)./ClustCell(1:2:end);
%%
print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\cluster - base sig');
%%
figure('name','AP cells in Cluster',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   420   420])

b =bar(1:10,sum(cluststim)/sum(sum(cluststim)));
b.LineStyle = 'none';
b.FaceColor = 'flat';
b.LineWidth = lnwd;
b.CData = piecol{1};

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.YLim = [0 1];
% ax.XLim = [0.5 1.5];
ax.LineWidth = lnwd;
% ax.YTick = 0:.2:1;
% ax.XTickLabel = {};

% title('AP cells in Cluster')
% ylabel('Clusters count')
ylabel('probability')
xlabel('AP cells per cluster')
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\clust AP cells')

%%
figure('name','Network Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   420   420])

b = bar(1:10,sum(clustplace)/sum(sum(clustplace)));
b.LineStyle = 'none';
b.FaceColor = 'flat';
b.LineWidth = lnwd;
b.CData = piecol{3};

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.YLim = [0 1];
% ax.XLim = [0.5 1.5];
ax.LineWidth = lnwd;
% ax.YTick = 0:.2:1;
% ax.XTickLabel = {};

% title('Place cells in cluster')
% ylabel('Clusters count')
ylabel('probability')
xlabel('PCs per cluster')
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\clust PCs')
load('cclustID.mat')
experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};
CAIMcorr = CAIM(:,:);

%%
figure('name','Cluster size',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   420   420])

b = bar(1:40,sum(clustsize)/sum(sum(clustsize)));
b.LineStyle = 'none';
b.FaceColor = 'flat';
b.LineWidth = lnwd;
b.CData = [1 1 1];

ax = gca;
ax.FontSize = ftsz;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.YLim = [0 1];
% ax.XLim = [0.5 1.5];
ax.LineWidth = lnwd;
% ax.YTick = 0:.2:1;
% ax.XTickLabel = {};

% title('Cluster size')
% ylabel('Clusters count')
ylabel('probability')
xlabel('Cells per cluster')
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\clust size')
%%
figure('name','cluster freq',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   420   420])

b = bar(.1:.1:1,sum(clustfreq)/sum(sum(clustfreq)));
b.LineStyle = 'none';
b.FaceColor = 'flat';
b.LineWidth = lnwd;
b.CData = [1 1 1];

ax = gca;
ax.FontSize = ftsz;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.YLim = [0 1];
% ax.XLim = [0.5 1.5];
ax.LineWidth = lnwd;
% ax.YTick = 0:.2:1;
% ax.XTickLabel = {};

% title('Cluster freq')
% ylabel('Clusters count')
ylabel('probability')
xlabel('% of net events')
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\clust freq')

%%
figure('name','Cell frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   220   620])
hold on

b = bar([100*a/sum(a); 100*a/sum(a)],'stacked');
for i = 1:length(b)
    b(i).LineStyle = 'none';
    b(i).FaceColor = 'flat';
    b(i).LineWidth = lnwd;
    b(i).CData = piecol{i};
end

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.YLim = [0 1];
ax.XLim = [0.5 1.5];
ax.LineWidth = lnwd;
% ax.YTick = 0:.2:1;
ax.XTickLabel = {};
ylabel('% of clusters')
% ylabel('participated net events %')
% xlabel('Ca^2^+events /min')

% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\cclust feature')

%% Clustersize AP vs PC vs None 

figure('name','Cluster size',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   420   420])
hold on
% b = bar(1:40,sum(clustsize)/sum(sum(clustsize)));
% b.LineStyle = 'none';
% b.FaceColor = 'flat';
% b.LineWidth = lnwd;
% b.CData = [1 1 1];
% plot(1:40,cumsum(sum(clustsize)/sum(sum(clustsize))));
plot(1:40,cumsum(sum(clustsizeAP)/sum(sum(clustsizeAP))),...
    'LineWidth',lnwd,...
    'color',piecol{1});
plot(1:40,cumsum(sum(clustsizePC)/sum(sum(clustsizePC))),...
    'LineWidth',lnwd,...
    'color',piecol{3});
plot(1:40,cumsum(sum(clustsizenone)/sum(sum(clustsizenone))),...
    'LineWidth',lnwd,...
    'color',piecol{4});
ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.YLim = [0 1];
% ax.XLim = [0.5 1.5];
ax.LineWidth = lnwd;
% ax.YTick = 0:.2:1;
% ax.XTickLabel = {};

% title('Cluster size')
% ylabel('Clusters count')
ylabel('fraction')
xlabel('Cells per cluster')

%% P value clusters

num = 1:8;
netp = nan(length(num),size(CAIM,2));
clustp = nan(length(num),size(CAIM,2));
netnum = nan(length(num),size(CAIM,2));
for i = num 
    
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
        if ~isempty(CAIM(i,k).network.clustsize)
            netp(i,k) = CAIM(i,k).network.netp(1);
            clustp(i,k) = CAIM(i,k).network.clustp;
            netnum(i,k) = length(CAIM(i,k).network.netpos);
        end
    end
    
end

subplot(1,2,1)
x = cumsum(histcounts(clustp(~isnan(clustp)),0:.05:1,'normalization','probability'));
x = [0 x];
plot(0:.05:1,x)
ylabel('cum prob')
xlabel('p-value')
grid on
subplot(1,2,2)
scatter(netnum(:),clustp(:))
grid on
ylabel('p value')
xlabel('net events in session')
%%
% subplot(1,3,3)
scatter(netp(:),clustp(:))
grid on
ylabel('p value cluster')
xlabel('p value network')
xlim([0 1])

% exportfig('C:\Users\Admin\Dropbox\matlab\clustp')

%% Shared place field and cluster

fieldclust = zeros(11,size(CAIM,2),6);
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
        kk = fullCAIM(j);      
        cclusttemp = CAIM(i,kk).cclust;
        isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
        isnet = cclusttemp(:,cclustID.netprob)>0;
        angles = cclusttemp(:,cclustID.plcvctang);
        clust = cclusttemp(:,cclustID.clust);
        
        % same field

        shiftangles = angles;
        shiftangles(angles<0) = shiftangles(angles<0)+2*pi;

        a = find(angles>0);
        samefield = zeros(length(angles));

        for jj = 1:length(a)
            for k = 1:length(a)
                if angles(a(jj))>=angles(a(k))-pi/4 && angles(a(jj))<=angles(a(k))+pi/4
                   samefield(a(jj),a(k)) = 1;          
                elseif shiftangles(a(jj))>=shiftangles(a(k))-pi/4 && shiftangles(a(jj))<=shiftangles(a(k))+pi/4
                   samefield(a(jj),a(k)) = 1;
                end
            end
        end
        samefield(1:length(samefield)+1:end) = 0;

        % same cluster

        sameclust = zeros(length(clust));
        b = find(clust);
        for jj = 1:length(b)
            for k = 1:length(b)
                if clust(b(jj)) == clust(b(k))
                   sameclust(b(jj),b(k)) = 1;          
                end
            end
        end

        sameclust(1:length(sameclust)+1:end) = 0;
        
        sameclustinvert = ones(size(sameclust))-sameclust;
        samefieldinvert = ones(size(sameclust))-samefield;
        samefieldinvert(1:length(samefield)+1:end) = 0;
        sameclustinvert(1:length(samefield)+1:end) = 0;
        %% read out values
%         isplace(a) =1;
        % 1 number of cells
        % 2 number of samefield cells
        % 3 number of samecluster cells
        % 4 number of place cells       
        % 5 number of shared cluster and direction
        % 6 number of shared cluster but not direction
        % 7 number of shared direction but not cluster
        % 8 number of shared nothing
        % 9 number of possible combinations
        % 10 number of all cells with vector
        % 11 number of all cells that share cluster and direction
        b = sameclust(isplace,isplace).*samefield(isplace,isplace) + sameclust(isplace,isplace).*samefieldinvert(isplace,isplace) + sameclustinvert(isplace,isplace).*samefield(isplace,isplace) +sameclustinvert(isplace,isplace).*samefieldinvert(isplace,isplace);
        
        fieldclust(:,kk,i-2) = [sum(~isnan(cclusttemp(:,cclustID.expID)));               
                sum(sum(samefield(isplace,isplace)))/2;
                sum(sum(sameclust(isplace,isplace)))/2;
                sum(isplace);              
                sum(sum(sameclust(isplace,isplace).*samefield(isplace,isplace)))/2;
                sum(sum(sameclust(isplace,isplace).*samefieldinvert(isplace,isplace)))/2;
                sum(sum(sameclustinvert(isplace,isplace).*samefield(isplace,isplace)))/2;
                sum(sum(sameclustinvert(isplace,isplace).*samefieldinvert(isplace,isplace)))/2;
                sum(sum(b))/2;
                length(a);
                sum(sum(sameclust(a,a).*samefield(a,a)))/2;];
        
%         figure
%         imagesc(b)
%             subplot(1,4,1)
%             imagesc(sameclust(isplace,isplace))
%             
%             subplot(1,4,2)
%             imagesc(sameclustinvert(isplace,isplace))
%             
%             subplot(1,4,3)
%             imagesc(samefield(isplace,isplace))
%             
%             subplot(1,4,4)
%             imagesc(samefieldinvert(isplace,isplace))
    end
end

%%
figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]],...
        'renderer','painters',...
        'visible','on')

x = 1:11;    
for i = 1:size(fieldclust,3)
    subplot(size(fieldclust,3),1,i)
    %%
%     i = 3;
    hold on
    pospair = .5*(fieldclust(4,:,i).*(fieldclust(4,:,i)-1));
%     pospair = ones(1,length(fieldclust(4,:,i)));

    a = fieldclust(5,:,i)./pospair;
    a(2,:) = fieldclust(6,:,i)./pospair;
    a(3,:) = fieldclust(7,:,i)./pospair;
    a(4,:) = fieldclust(8,:,i)./pospair;
%     a(5,:) = fieldclust(9,:,i)./pospair;
%     a(:,9) = NaN;
%     a(:,end+1) = nanmean(a,2); 
    a(:,end+1) = nansum(fieldclust(5:8,:,i),2)/sum(pospair);
    
    a(isnan(a)) = 0;
    a(a==inf | a==-inf) = 0;
    
    bar(x,a','stacked')
    if i == 1
        legend({'field & clust' 'clust' 'field' 'none'},'Location','best','Orientation','horizontal')
%         legend('boxoff')
    end
    
    ax = gca;
    ax.XTick = x;
%     ax.XLim = [0 12];
    ax.XTickLabel = [mouse(:); 'pool'];
    title(experiment(i+2))
    hold off
end