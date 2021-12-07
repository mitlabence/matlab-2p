% addpath(genpath('C:\Users\Admin\Documents\MATLAB\sorter'))
% addpath(genpath('C:\Users\Admin\Documents\MATLAB\treadmill'))
addpath(genpath('C:\Users\Admin\Documents\MATLAB\ca_source_extraction'))%
% Cues = load('Z:\Martin Pofahl\Cues.mat');
Baseline = load('\');
% Airpuff = load('Z:\Martin Pofahl\Airpuff.mat');
%load('Z:\Martin Pofahl\BigFatCluster.mat');
%load('cclustID.mat')
%%
lnwd = 3;
ftsz = 25;

%% Number of rounds
figure('name','number of rounds',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])

hold on
mouseID = Baseline.Behave.mouseID;
yin = Baseline.Behave.numrounds;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

% y2 = Cues.Fire.placecode(cell2mat(Cues.Fire.mouseID(:,2))==trial,j);
mouseID = Cues.Behave.mouseID;
yin = Cues.Behave.numrounds;
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Behave.mouseID;
yin =  Airpuff.Behave.numrounds;
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y y1 y2];
[~,~,stats] = anova1(ystat,'','off');
c = multcompare(stats,'Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];

x = [1 3 5];
b = bar(x, y);
b.FaceColor = 'flat';
% b.LineWidth = lnwd;
b.CData(1,:) = [0 176/255,240/255];
b.CData(2,:) = [0 176/255,240/255];
b.CData(3,:) = [0 176/255,240/255];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',2)

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 6];
ax.LineWidth =lnwd;
ax.XTick = [1 3 5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};

ylabel('# rounds')

hold off 
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\rounds')

%% Runtime 
figure('name','Plce field',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])

hold on
mouseID = Baseline.Behave.mouseID;
yin = Baseline.Behave.runtime/60;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

% y2 = Cues.Fire.placecode(cell2mat(Cues.Fire.mouseID(:,2))==trial,j);
mouseID = Cues.Behave.mouseID;
yin = Cues.Behave.runtime/60;
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Behave.mouseID;
yin =  Airpuff.Behave.runtime/60;
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y y1 y2];
[~,~,stats] = anova1(ystat,'','off');
c = multcompare(stats,'Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];

x = [1 3 5];
b = bar(x, y);
b.FaceColor = 'flat';
% b.LineWidth = lnwd;
b.CData(1,:) = [.4 176/255,120/255];
b.CData(2,:) = [.4 176/255,120/255];
b.CData(3,:) = [.4 176/255,120/255];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',2)

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 6];
ax.LineWidth =lnwd;
ax.XTick = [1 3 5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};

ylabel('Running time / min')

hold off 
exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\runtime')

%% pupil size bar graph

figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Behave.mouseID;
yin = Baseline.Behave.pupilsize(1,2:3,:);
yin = yin./Baseline.Behave.pupilsize(3,2:3,:);
yin = permute(yin,[3,2,1]);
yin(cell2mat(mouseID(:,2))==1,:) = NaN;
yin(cell2mat(mouseID(:,2))==2,:) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end
% y1 = Cues.Fire.meanfire(cell2mat(Cues.Fire.mouseID(:,2))==trial,2:3);
mouseID = Cues.Behave.mouseID;
yin = Cues.Behave.pupilsize(1,2:3,:);
yin = yin./Cues.Behave.pupilsize(3,2:3,:);
yin = permute(yin,[3,2,1]);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Behave.mouseID;
yin =  Airpuff.Behave.pupilsize(1,2:3,:);
yin =  yin./Airpuff.Behave.pupilsize(3,2:3,:);
yin = permute(yin,[3,2,1]);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y(~isnan(y(:,1)),:); y1(~isnan(y(:,1)),:); y2(~isnan(y(:,1)),:)];
reps = size(y(~isnan(y(:,1)),:),1);
[~,~,stats] = anova2(ystat,reps,'off');
c = multcompare(stats,'Display','off');
c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];

b = bar(x, y);
b.LineStyle = 'none';
b.FaceColor = 'flat';
b.LineWidth = lnwd;
b.CData(1,:) = [1 1 0];
b.CData(2,:) = [.5 .5 0];
b.CData(3,:) = [1 1 0];
b.CData(4,:) = [.5 .5 0];
b.CData(5,:) = [1 1 0];
b.CData(6,:) = [.5 .5 0];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',lnwd);

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.YLim = [.8 1];
ax.LineWidth = lnwd;
ax.XLim = [0 9];
ax.XTick = [1.5 4.5 7.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
% ax.XTickLabel = {};
ylabel('Pupil size/max','FontSize',ftsz-2)
        
% [h,p] = ttest(Fire.meanfire(:,2),Fire.meanfire(:,3));
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\pupil size')

%% bargraph mean firing frequency of active cells
figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Fire.mouseID;
yin = Baseline.Fire.meanfire(:,2:3);
yin(cell2mat(mouseID(:,2))==1,:) = NaN;
yin(cell2mat(mouseID(:,2))==2,:) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end
% y1 = Cues.Fire.meanfire(cell2mat(Cues.Fire.mouseID(:,2))==trial,2:3);
mouseID = Cues.Fire.mouseID;
yin = Cues.Fire.meanfire(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Fire.mouseID;
yin =  Airpuff.Fire.meanfire(:,2:3);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y; y1 ; y2];
reps = size(y,1);
[~,~,stats] = anova2(ystat,reps,'off');
c = multcompare(stats,'Display','off');
c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];

y = [nanmean(y) nanmean(y1) nanmean(y2)];

b = bar(x, y);
b.LineStyle = 'none';
b.FaceColor = 'flat';
b.LineWidth = lnwd;
b.CData(1,:) = [.6 .6 .6];
b.CData(2,:) = [.2 .2 .2];
b.CData(3,:) = [.6 .6 .6];
b.CData(4,:) = [.2 .2 .2];
b.CData(5,:) = [.6 .6 .6];
b.CData(6,:) = [.2 .2 .2];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',lnwd);

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 6];
ax.YLim = [0 2.5];
ax.LineWidth = lnwd;
ax.XLim = [0 9];
ax.XTick = [1.5 4.5 7.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
% ax.XTickLabel = {};
ylabel('Events min^-^1','FontSize',ftsz-2)
        
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\meanfire box')

%% bargraph mean amplitude of active cells
figure('name','Firing Amplitude',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Fire.mouseID;
jj = 1; %1 for mean, 2 for std, 3 for max
j = [ 2 3];
yin = Baseline.Fire.meanamplitude(jj:3:end,j);
yin(cell2mat(mouseID(:,2))==1,:) = NaN;
yin(cell2mat(mouseID(:,2))==2,:) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end
% y1 = Cues.Fire.meanfire(cell2mat(Cues.Fire.mouseID(:,2))==trial,2:3);
mouseID = Cues.Fire.mouseID;
yin = Cues.Fire.meanamplitude(jj:3:end,j);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Cues.Fire.mouseID;
yin = Airpuff.Fire.meanamplitude(jj:3:end,j);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

ystat = [y; y1; y2];
reps = size(y,1);
[~,~,stats] = anova2(ystat,reps,'off');
c = multcompare(stats,'Display','off');
c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

yerr = [std(y)/sqrt(length(y)) std(y1)/sqrt(length(y1)) std(y2)/sqrt(length(y2))];

y = [mean(y) mean(y1) mean(y2)];

b = bar(x, y);
b.LineStyle = 'none';
b.FaceColor = 'flat';
b.LineWidth = lnwd;
b.CData(1,:) = [.6 .6 .6];
b.CData(2,:) = [.2 .2 .2];
b.CData(3,:) = [.6 .6 .6];
b.CData(4,:) = [.2 .2 .2];
b.CData(5,:) = [.6 .6 .6];
b.CData(6,:) = [.2 .2 .2];
errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',lnwd);

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 9];
ax.YLim = [0 .8];
ax.LineWidth = lnwd;
ax.XTick = [1.5 4.5 7.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
ylabel('\DeltaF/F','FontSize',ftsz-2) 

% [h,p] = ttest(Fire.meanamplitude(jj:3:end,j(1)),Fire.meanamplitude(jj:3:end,j(2)));
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\meanamp')

%% Bar graph number of network events

figure('name','Network Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Fire.mouseID;
yin = Baseline.Network.netfreq(:,2:3);
yin(cell2mat(mouseID(:,2))==1,:) = NaN;
yin(cell2mat(mouseID(:,2))==2,:) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end
% y1 = Cues.Fire.meanfire(cell2mat(Cues.Fire.mouseID(:,2))==trial,2:3);
mouseID = Cues.Fire.mouseID;
yin =  Cues.Network.netfreq(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Fire.mouseID;
yin =  Airpuff.Network.netfreq(:,2:3);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y; y1; y2];
reps = size(y,1);
[~,~,stats] = anova2(ystat,reps,'off');
c = multcompare(stats,'Display','off');
c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];
b = bar(x, y);
hold on
b.FaceColor = 'flat';
b.LineStyle = 'none';
b.CData(1,:) = [.1 1 .1];
b.CData(2,:) = [0 .5 0];
b.CData(3,:) = [.1 1 .1];
b.CData(4,:) = [0 .5 0];
b.CData(5,:) = [.1 1 .1];
b.CData(6,:) = [0 .5 0];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',lnwd)

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 9];
ax.YLim = [0 ceil(max(y)+max(yerr))];
ax.LineWidth = lnwd;
ax.XTick = [1.5 4.5 7.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
% ax.XTickLabel = {};
ax.YTick = 0:2:6;
ylabel('Net min^-^1','FontSize',ftsz-2)
box('off')
hold off

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\netfreq')
%% Bar graph of network occurence

figure('name','Network Occurence',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Fire.mouseID;
yin = Baseline.Network.netoccur(:,2:3);
yin(cell2mat(mouseID(:,2))==1,:) = NaN;
yin(cell2mat(mouseID(:,2))==2,:) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end
% y1 = Cues.Fire.meanfire(cell2mat(Cues.Fire.mouseID(:,2))==trial,2:3);
mouseID = Cues.Fire.mouseID;
yin =  Cues.Network.netoccur(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Fire.mouseID;
yin =  Airpuff.Network.netoccur(:,2:3);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y; y1; y2];
reps = size(y,1);
[~,~,stats] = anova2(ystat,reps,'off');
c = multcompare(stats,'Display','off');
c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];
b = bar(x, y);
hold on
b.FaceColor = 'flat';
b.LineStyle = 'none';
b.CData(1,:) = [.1 1 .1];
b.CData(2,:) = [0 .5 0];
b.CData(3,:) = [.1 1 .1];
b.CData(4,:) = [0 .5 0];
b.CData(5,:) = [.1 1 .1];
b.CData(6,:) = [0 .5 0];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',lnwd)

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 9];
ax.YLim = [0 ceil(max(y)+max(yerr))];
ax.LineWidth = lnwd;
ax.XTick = [1.5 4.5 7.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
ylabel('Net session^-^1','FontSize',ftsz-2)
box('off')
hold off

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\netoccur')

%% Network size cumulative
figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   350   420])

a1 = [];
a2 = [];
a3 = [];
for i = 1:size(CAIM,2)
    %Baseline   
    for j = 3:5
        if ~isempty(CAIM(j,i).A)
            a1 = [a1 CAIM(j,i).network.netnum];
        end
    end
    %Cues   
    for j = 6:8
        if ~isempty(CAIM(j,i).A)
            a2 = [a2 CAIM(j,i).network.netnum];
        end
    end
    
    %AP
    for j = 9:11
        if ~isempty(CAIM(j,i).A)
            a3 = [a3 CAIM(j,i).network.netnum];
        end
    end
end

ystat = [a1 a2 a3];
ygroup(1:length(a1)) = 1;
ygroup(length(a1)+1:length(a1)+length(a2)) = 2;
ygroup(length(a1)+length(a2)+1:length(a1)+length(a2)+length(a3)) = 3;
[~,~,stats] = kruskalwallis(ystat,ygroup,'off');
c = multcompare(stats,'Display','off');

[aa1,bb1]  = histcounts(a1, 'Normalization', 'probability'); 
[aa2,bb2]  = histcounts(a2, 'Normalization', 'probability');
[aa3,bb3]  = histcounts(a3, 'Normalization', 'probability');
          
plot(bb1(1:end-1)+.5,cumsum(aa1),...
    'linewidth',lnwd,...
    'color',[.7 .2 .5])
hold on
plot(bb2(1:end-1)+.5,cumsum(aa2),...
    'linewidth',lnwd,...
    'color',[.3 .4 .8])
plot(bb3(1:end-1)+.5,cumsum(aa3),...
    'linewidth',lnwd,...
    'color',[.2 .8 .5])

ax = gca;
ax.XLim = [2.5 30]; 
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.LineWidth = lnwd;
box off
xlabel('Cells event^-^1','FontSize',ftsz-2)
ylabel('probability','FontSize',ftsz-2)
hold off

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\netsize')

%% Trigered network plot
figure('name','NEtwork events',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   350   420])

stim = Baseline.Response.network;
% stim = Response.network;
win = 01:121;

hold on
sumresp = zeros(size(stim.sumresp));
for i = 1:size(stim.nums,1)
    sumresp(i,:) = stim.sumresp(i,:)/stim.nums(i,1);%/stim.nums(i,2);
end

% plot([0 0],[0 .31],'g','linewidth',2)
a = mean(stim.bulk(:,win));
b = std(stim.bulk(:,win))/sqrt(size(stim.bulk(:,win),1));

bar(stim.times(1,win)/1000,mean(sumresp(:,win)),'g',...    'EdgeColor','g',...
    'linewidth',.1*lnwd);

fill([stim.times(1,win)/1000,fliplr(stim.times(1,win)/1000)],[2*a-2*b,fliplr(2*a+2*b)]+1.5,[1 0 0],...
    'EdgeColor',[1 0 0],...
    'EdgeAlpha',0,...
    'FaceAlpha',.5)

plot(stim.times(1,win)/1000,2*a+1.5,'r',...
    'linewidth',lnwd-.5)

ax = gca;
% ax.FontSize = ftsz-2;
% ax.YColor = [1 1 1];
% ax.XColor = [1 1 1];
% ax.Color = [0 0 0];
ax.XLim = [-1 2];
% ax.YLim = [0 1];
% ax.LineWidth = lnwd;

% xlabel('time s^-^1','FontSize',ftsz-2)
% ylabel('probability','FontSize',ftsz-2)

% axis off
hold off

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\MPP Net')

%% Bar graph MPP amplitude
figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Bulk.mouseID;
yin = Baseline.Bulk.bulkbase(:,2:3);
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Cues.Bulk.mouseID;
yin = Cues.Bulk.bulkbase(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Bulk.mouseID;
yin =  Airpuff.Bulk.bulkbase(:,2:3);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y; y1; y2];
reps = size(y,1);
[~,~,stats] = anova2(ystat,reps,'off');
c = multcompare(stats,'Display','off');
c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];

b = bar(x, y);
b.FaceColor = 'flat';
b.LineStyle = 'none';
b.LineWidth = lnwd ;
b.CData(1,:) = [1 .2 .2];
b.CData(2,:) = [.5 0 0];
b.CData(3,:) = [1 .2 .2];
b.CData(4,:) = [.5 0 0];
b.CData(5,:) = [1 .2 .2];
b.CData(6,:) = [.5 0 0];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',lnwd)

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [ 1 1 1];
ax.Color = [0 0 0];
ax.YLim = [-.5 1.1];
ax.LineWidth = lnwd;
ax.XLim = [0 9];
ax.XTick = [1.5 4.5 7.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
% ax.XTickLabel = {};
ylabel('Z score','FontSize',ftsz-2) 

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\MPP Ampl')

%% Bar graph MPP standart deviation
figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

% Average over trials
mouseID = Baseline.Bulk.mouseID;
yin = Baseline.Bulk.bulkstd(:,2:3);
yin([1 12 14],2) = 0; % correct for non-runner mice
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & yin(:,2)>0,:),1);
end

mouseID = Cues.Bulk.mouseID;
yin = Cues.Bulk.bulkstd(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Bulk.mouseID;
yin =  Airpuff.Bulk.bulkstd(:,2:3);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y; y1; y2];
reps = size(y,1);
[~,~,stats] = anova2(ystat,reps,'off');
c = multcompare(stats,'Display','off');
c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];

b = bar(x, y);
b.FaceColor = 'flat';
b.LineStyle = 'none';
b.CData(1,:) = [1 .2 .2];
b.CData(2,:) = [.5 0 0];
b.CData(3,:) = [1 .2 .2];
b.CData(4,:) = [.5 0 0];
b.CData(5,:) = [1 .2 .2];
b.CData(6,:) = [.5 0 0];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',lnwd)

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 9];
ax.YLim = [0 1.2];
ax.LineWidth = lnwd;
ax.XTick = [1.5 4.5 7.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
ylabel('\sigma','FontSize',ftsz-2)        

hold off 

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\MPP std')

%% Correlation MPP and speed

num = 3:8;
speedcorr = nan(6,10,3);
for i = 1:length(num)%1:size(CAIM,1)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).bulk)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    
    for j = 1:length(fullCAIM)
        k = fullCAIM(j);    
        speedcorr(i,k,:) = CAIM(num(i),k).bulk.speedcorr;

    end
end

pSpeed = speedcorr(:,:,2);
rSpeed = speedcorr(:,:,1);
rSpeed = rSpeed(:,~isnan(rSpeed(1,:)));
rSpeed = [nanmean(rSpeed(1:3,:),1); nanmean(rSpeed(6,:),1)];
figure
plot(rSpeed)
hold on
boxplot(rSpeed')
% anova1(rSpeed')


%% Correlation of MPP and network

figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   350   420])
hold on

% Average over trials
mouseID = Baseline.Bulk.mouseID;
yin = Baseline.Bulk.netxbulk(2:2:end,:);
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

a = mean(y);
b = std(y)./sqrt(size(y,1));

fill([Baseline.Bulk.netxbulk(1,:),fliplr(Baseline.Bulk.netxbulk(1,:))],[a-b,fliplr(a+b)],[.7 .7 0],...
    'EdgeColor',[.7 .7 .7],...
    'EdgeAlpha',1,...
    'FaceAlpha',.5)
hold on

plot([Baseline.Bulk.netxbulk(1,1) Baseline.Bulk.netxbulk(1,end)],[0 0],':','color',[.5 .5 .5])
plot([0 0],[-2000 4000],':','color',[.5 .5 .5])
plot(Baseline.Bulk.netxbulk(1,:),a,...
    'Color',[1 1 0],...
    'LineWidth',lnwd)
hold off
box off

ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [-4 4];
% ax.XTick = -5:5;
% ax.XTickLabel = {};
ax.YLim = [-1000 3000];
ax.YTick = ax.YLim(1):1000:ax.YLim(2);
ax.YTickLabel = ax.YTick/1000;
ax.LineWidth = lnwd;
xlabel('\Deltat s^-^1','linewidth',lnwd)

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\MPP corr')

%% Number of place fields
figure('name','Place field',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   450   420])

hold on
j=2;
mouseID = Baseline.Fire.mouseID;
yin = Baseline.Fire.placecode(:,j)./Baseline.Fire.numcells;
yin(yin==0) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

% y2 = Cues.Fire.placecode(cell2mat(Cues.Fire.mouseID(:,2))==trial,j);
mouseID = Cues.Fire.mouseID;
yin = Cues.Fire.placecode(:,j)./Cues.Fire.numcells;
yin(yin==0) = NaN;
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Fire.mouseID;
yin =  Airpuff.Fire.placecode(:,j)./Airpuff.Fire.numcells;
yin(yin==0) = NaN;
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y y1 y2];
[~,~,stats] = anova1(ystat,'','off');
c = multcompare(stats,'Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];

x = [1 3 5];
b = bar(x, y);
b.FaceColor = 'flat';
% b.LineWidth = lnwd;
b.CData(1,:) = [0 176/255,240/255];
b.CData(2,:) = [0 176/255,240/255];
b.CData(3,:) = [0 176/255,240/255];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',2)

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 6];
ax.LineWidth =lnwd;
ax.XTick = [1 3 5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};

ylabel('fraction of PCs %')

hold off 
exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\place field fraction')

%% Size of place fields
figure('name','Plce field size',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   450   420])

hold on
% 1 mean of vector, 2 # significant place vectors, 3 # of significant place fields 
j=4;
mouseID = Baseline.Fire.mouseID;
yin = Baseline.Fire.placecode(:,j);
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

% y2 = Cues.Fire.placecode(cell2mat(Cues.Fire.mouseID(:,2))==trial,j);
mouseID = Cues.Fire.mouseID;
yin = Cues.Fire.placecode(:,j);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Fire.mouseID;
yin =  Airpuff.Fire.placecode(:,j);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y y1 y2];
[~,~,stats] = anova1(ystat,'','off');
c = multcompare(stats,'Display','off');

yerr = [nanstd(y)/sqrt(length(y)) nanstd(y1)/sqrt(length(y1))  nanstd(y2)/sqrt(length(y2))];
y = [nanmean(y) nanmean(y1) nanmean(y2)];

x = [1 3 5];
b = bar(x, y);

b.FaceColor = 'flat';
% b.LineWidth = lnwd;
b.CData(1,:) = [0 176/255,240/255];
b.CData(2,:) = [0 176/255,240/255];
b.CData(3,:) = [0 176/255,240/255];

errorbar(x,y,yerr,'.',...
            'Color',[1 1 1],...
            'LineWidth',2)

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XLim = [0 6];
ax.LineWidth =lnwd;
ax.XTick = [1 3 5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};

ylabel('place field size /cm')   

hold off 
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\place field size')


%% place fields examples
figure('name','Plce fields',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   450   720])

use = [5 8]; % exp trials to compare
mice = [103 155 158 194 195 224 226 227];% Mice to use cells from

cellID = cclust(:,cclustID.expID,:);
cellID = max(cellID,[],3);
cellID = cellID == mice;
cellID = max(cellID,[],2);
thresh = 1;% threshold of appearance

samecelltemp = samecell(use,cellID);
actcell = samecelltemp;
actcell(actcell>0) = 1;
actcell = sum(actcell,1);

cclusttemp = cclust;
exclude = samecell' == 0;
for i = 1:size(cclusttemp,3)
    cclusttemp(exclude(:,i),:,i) = NaN;
end
cclusttemp = cclust(cellID,:,:);
cclusttemp = cclusttemp(actcell == 1,:,use);
samecelltemp = samecelltemp(:,actcell == 1);



a = cclusttemp(:,cclustID.plcvct,:)>0 & cclusttemp(:,cclustID.plcvctp,:)<=.05;
b = cclusttemp(max(a(:,1,:),[],3),:,:);

plcfieldtemp = plcfield(cellID,:,:);
plcfieldtemp = plcfieldtemp(actcell == 1,:,use);
plcfieldtemp = plcfieldtemp(max(a(:,1,:),[],3),:,:); 
for i = 1:size(plcfieldtemp,1)
    for j = 1:size(plcfieldtemp,3)
        plcfieldtemp(i,:,j) = (plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)))./max(plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)));
    end
end

% plcfields in baseline condition
% sort angles
plccenter = b(:,cclustID.plcvctang,1);%plcfld
plccenter = plccenter(b(:,cclustID.plcvct,1)>0);
[plccenter,aa] = sort(plccenter);
plcfield1 = plcfieldtemp(b(:,cclustID.plcvct,1)>0,:,1);
plcfield1 = plcfield1(aa,:);
c = find(abs(plccenter)==min(abs(plccenter)))-7;
plcfield1 = plcfield1([c:end 1:c-1],:);

subplot(2,1,1)
imagesc(plcfield1)

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.XLim = [0 4];
ax.LineWidth = lnwd;
ax.XTick = [];
ax.XTickLabel = [];
ylabel('cell ID')

% plcfields in basline condition
% sort angles
plccenter = b(:,cclustID.plcvctang,2);%plcfld
plccenter = plccenter(b(:,cclustID.plcvct,2)>0);
[plccenter,aa] = sort(plccenter);
plcfield2 = plcfieldtemp(b(:,cclustID.plcvct,2)>0,:,2);
plcfield2 = plcfield2(aa,:);
c = find(abs(plccenter)==min(abs(plccenter)))-6;
plcfield2 = plcfield2([c:end 1:c-1],:);


subplot(2,1,2)
imagesc(plcfield2)

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.XLim = [0 4];
ax.LineWidth = lnwd;
ax.XTick = [33 66 100];
ax.XTickLabel = [50:50:150];
ylabel('cell ID')
xlabel('Belt /cm')

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\place fields')

%% pie chart of cell feature

trial = 11;
cclusttemp = [];
for i = 1:size(CAIM,2)
    for j = 1:length(trial)
        k = trial(j);
        cclusttemp = [cclusttemp; CAIM(k,i).cclust];
    end
end

isnet = cclusttemp(:,cclustID.netprob)>0;
isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
isplace = false(size(isplace));
isstim = cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05;
isstim = false(size(isstim));

a = zeros(1,8);
% network, no place, stim
a(1) = sum(isnet &  ~isplace & isstim)/length(cclusttemp);
% network, place, stim
a(2) = sum(isnet &  isplace & isstim)/length(cclusttemp);
% network, place, no stim
a(3) = sum(isnet &  isplace & ~isstim)/length(cclusttemp);
% network, no place, no stim
a(4) = sum(isnet  &  ~isplace  & ~isstim)/length(cclusttemp);
% no network, no place, stim
a(5) = sum(~isnet &  ~isplace  & isstim)/length(cclusttemp);
% no network, place, stim
a(6) = sum(~isnet &  isplace & isstim)/length(cclusttemp);
% no network, place, no stim
a(7) = sum(~isnet &  isplace & ~isstim)/length(cclusttemp);
% no network, no place, no stim
a(8) = sum(~isnet &  ~isplace  & ~isstim)/length(cclusttemp);
a = (100*a);

% labels = {'net & stim','net & plc & stim','net & plc','net','stim','plc & stim','plc','non'};
labels = {'','','','','','','',''};
piecol = {[1 0 1] , [1 1 0]   ,[0 1 1] ,[0 1 0]  ,[.5 0 .5]        ,[.5 .5 0]      ,[0 .5 .5],[.5 .5 .5]};
labels = labels(a>0);
piecol = piecol(a>0);
a = a(a>0);
if ~isempty(find(a,1))
    figure('name',['Feature pie ' num2str(trial)],...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   450   450])

    p = pie(a,labels);
    for i = 1:length(a)
        pp = p(2*i-1);
        pp.FaceColor  = piecol{i};
    end
end
size(cclusttemp,1)
% exportfig(['C:\Users\Admin\Dropbox\Lab\Vortr�ge\feature pie ' num2str(trial)])
%% number of events of active cells (histogram)
figure('name','Cell frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   350   420])

hold on
trial = 9:11;
cclusttemp = [];
for i = 1:size(CAIM,2)
    for j = 1:length(trial)
        k = trial(j);
        cclusttemp = [cclusttemp; CAIM(k,i).cclust];
    end
end

isnet = cclusttemp(:,cclustID.netprob)>0;
isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
% isplace = false(size(isplace));
isstim = cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05;

fplace = cclusttemp(isplace,cclustID.meanf);
fstim = cclusttemp(isstim,cclustID.meanf);
fother = cclusttemp(~(isplace|isstim),cclustID.meanf);

maxcount = 50;
%xtil = string(0:5:maxcount);
%xtil(end) = string([' >' num2str(maxcount)]); 
%fireprob1(fireprob1>maxcount) = maxcount;
[aa1,bb1]  = histcounts(fother,...
    'BinWidth',.5,...
    'BinLimits',[-1, maxcount],...
    'Normalization', 'probability'); 
plot(bb1(1:end-1)+.5,cumsum(aa1),...
    'linewidth',lnwd,...
    'color',[0 1 0])
[aa2,bb2]  = histcounts(fplace,...
    'BinWidth',.5,...
    'BinLimits',[-1, maxcount],...
    'Normalization', 'probability');
plot(bb2(1:end-1)+.5,cumsum(aa2),...
    'linewidth',lnwd,...
    'color',[0 1 1])
[aa3,bb3]  = histcounts(fstim,...
    'BinWidth',1,...
    'BinLimits',[-1, maxcount],...
    'Normalization', 'probability');
plot(bb3(1:end-1)+.5,cumsum(aa3),...
    'linewidth',lnwd,...
    'color',[1 0 1])

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
% ax.XLim = [0 4];
ax.LineWidth = lnwd;
% ax.XTick = [33 66 100];
% ax.XTickLabel = [50:50:150];
ylabel('probability')
xlabel('Ca^2^+events /min')
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\events cells')

%% Correlation netprob and meanf
figure('name','Cell frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   420   420])
hold on

a = cclusttemp(isnet,cclustID.meanf);
c = cclusttemp(isnet & isstim,cclustID.meanf);
d = cclusttemp(isnet & isplace,cclustID.meanf);
e = cclusttemp(isnet & isplace & isstim,cclustID.meanf);

b = cclusttemp(isnet,cclustID.netprob);
cc = cclusttemp(isnet & isstim,cclustID.netprob);
ee = cclusttemp(isnet & isplace & isstim,cclustID.netprob);
dd = cclusttemp(isnet & isplace,cclustID.netprob);

% b = cclusttemp(isnet,cclustID.netperc);
% cc = cclusttemp(isnet & isstim,cclustID.netperc);
% ee = cclusttemp(isnet & isplace & isstim,cclustID.netperc);
% dd = cclusttemp(isnet & isplace,cclustID.netperc);

scatter(a,b,'*','MarkerEdgeColor',[0 1 0])
scatter(c,cc,'*','MarkerEdgeColor',[1 0 1])
scatter(d,dd,'*','MarkerEdgeColor',[0 1 1])
scatter(e,ee,'*','MarkerEdgeColor',[1 1 0])

bin = 100;
binb = zeros(bin-1,1);
bina = zeros(bin-1,1);
errora = zeros(bin-1,1);
errorb = zeros(bin-1,1);
bin = (min(a):(max(a)-min(a))/bin:max(a))';
for i = 1:length(bin)-1
    bina(i) = nanmean(a(a>=bin(i) & a<=bin(i+1)));
    binb(i) = nanmean(b(a>=bin(i) & a<=bin(i+1)));
    errora(i) = nanstd(a(a>=bin(i) & a<=bin(i+1)));
    errorb(i) = nanstd(b(a>=bin(i) & a<=bin(i+1)))/sqrt(length(b(a>=bin(i) & a<=bin(i+1))));
end
% errorbar(bina,binb,errorb/2,errorb/2,errora/2,errora/2,'*',...
%     'Marker','none',...
%     'color',[1 0 0],...
%     'LineWidth',lnwd-1)

[h,p] = corr(a,b);
% title(['Network against activity, r = ' num2str(h) ', p = ' num2str(p)]);

ax = gca;
ax.FontSize = ftsz-5;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.YLim = [0 1];
ax.XLim = [0 50];
ax.LineWidth = lnwd;
ax.YTick = 0:.2:1;
ax.YTickLabel = 0:20:100;
ylabel('network related events %')
% ylabel('participated net events %')
xlabel('Ca^2^+events /min')

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\corr 2')

%% AP response pool plot
figure('name','Cell frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   520   620])

stim = Airpuff.Response.airpuff;
win = 16:91;
sumresp = zeros(size(stim.sumresp));
nresp = zeros(size(stim.nresp));
trial = 1;
for i = 1:size(stim.nums,1)
    sumresp(i,:) = stim.sumresp(i,:)/stim.nums(i,1);%/stim.nums(i,2);
    nresp(i,:) = stim.nresp(i,:)/stim.nums(i,1);%/stim.nums(i,2);
end
% sumresp = sumresp(cell2mat(stim.mouseID(:,2))==3,:);

exclude = [ ];
for i = 1:length(exclude)
    temp = find(cell2mat(stim.mouseID(:,4))==exclude(i));
    for j = 1:length(temp)
        stim.mouseID{temp(j),2} = 0;
    end
end
excludebulk = [ ];
for i = 1:length(excludebulk)
    temp = find(cell2mat(stim.mouseIDbulk(:,4))==excludebulk(i));
    for j = 1:length(temp)
        stim.mouseIDbulk{temp(j),2} = 0;
    end
end

x = stim.times(1,win)/1000;

subplot(5,1,1)

y = mean(sumresp(cell2mat(stim.mouseID(:,2))==trial,win),1);
yy = mean(nresp(cell2mat(stim.mouseID(:,2))==trial,win),1);
plot(x,y,'color',[0 1 0],'linewidth',lnwd-1)
hold on
b = std(sumresp(cell2mat(stim.mouseID(:,2))==trial,win),1)./sqrt(size(sumresp(cell2mat(stim.mouseID(:,2))==trial,win),1));

fill([x,fliplr(x)],[y-b,fliplr(y+b)],[0 .7 0],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)

plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd-1)
axis tight
axis off
hold off

subplot(5,1,2)
plot(x,yy,'color',[0 .7 0],'linewidth',lnwd-1)
hold on

b = std(nresp(cell2mat(stim.mouseID(:,2))==trial,win),1)./sqrt(size(nresp(cell2mat(stim.mouseID(:,2))==trial,win),1));
fill([x,fliplr(x)],[yy-b,fliplr(yy+b)],[0 .7 0],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)

plot([0 0],[min(yy-b) max(yy+b)],'m','linewidth',lnwd-1)
axis tight
axis off
hold off

subplot(5,1,3)
y = mean(stim.bulk(cell2mat(stim.mouseIDbulk(:,2))==trial,win),1);
b = std(stim.bulk(cell2mat(stim.mouseIDbulk(:,2))==trial,win),1)./sqrt(size(stim.bulk(cell2mat(stim.mouseIDbulk(:,2))==trial,win),1));
plot(x,y,'r','linewidth',lnwd-1)
hold on
fill([x,fliplr(x)],[y-b,fliplr(y+b)],[.7 0 0],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)
plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd-1)
axis tight
axis off
hold off

subplot(5,1,4)
y = stim.speed(cell2mat(stim.mouseID(:,2))==trial,win);
b = std(y,1)./sqrt(size(y,1));
y = mean(y,1);
plot(x,y,'color',[0 176/255,240/255],'linewidth',lnwd-1)
hold on
fill([x,fliplr(x)],[y-b,fliplr(y+b)],[0 176/255,240/255],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)
plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd-1)
axis tight
axis off
hold off

subplot(5,1,5)
y = stim.pupil(cell2mat(stim.mouseID(:,2))==trial,win);
b = std(y,1)./sqrt(size(y,1));
y = mean(y,1);
plot(x,y,'y','linewidth',lnwd-1)
hold on
fill([x,fliplr(x)],[y-b,fliplr(y+b)],'y',...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.2)
plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd-1)
axis tight
axis off
hold off

exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\AP resp pool')


%% cumsum of net number for complete data set
experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};

meannet = NaN(size(CAIM));
numclust = NaN(size(CAIM));
numcell = NaN(size(CAIM));
netp = NaN(size(CAIM));

figure('color',[1 1 1],...
        'position',[500 50 2*[594 420]],...
        'renderer','painters',...
        'visible','on')
        
for j = 1:size(CAIM,1) 
    subplot(4,4,j)
    for i = 1:size(CAIM,2)  
        if ~isempty(CAIM(j,i).A) && size(CAIM(j,i).network.netraster,2)>5
            [aa,bb] = histcounts((CAIM(j,i).network.netnum), 'Normalization', 'probability');
            numcell(j,i) = size(CAIM(j,i).network.netraster,1);
%             plot(bb(1:end-1),aa,'color',[numcell(j,i)/600 0 1-numcell(j,i)/600])
            
            plot(bb(1:end-1)+.5,cumsum(aa),'color',[numcell(j,i)/590 0 1-numcell(j,i)/590])
            hold on
            sumaa = 0;
%             meannet(j,i) = bb(find(aa==max(aa),1))+.5;
            for ii = 1:length(aa)
                sumaa = sumaa + aa(ii);
                if sumaa >0.66
                    meannet(j,i) = bb(ii)+.5;
                    break
                end
            end            
            numclust(j,i) = round(numcell(j,i)/meannet(j,i));
%             netp(j,i) = CAIM(j,i).network.netp(1);
        end
    end
    ax = gca;
    ax.XLim = [3 40];  
    ax.YLim = [0 1]; 
    title(experiment{j});
    xlabel('Cells per network');
    ylabel('fract of networks')
end

%printpdf('Network size')

%% Network median vs number off cells
% hold on
figure('color',[0 0 0],...
        'position',[500 100 2*[594 220]],...
        'renderer','painters',...
        'visible','on')
    
subplot(1,3,1)
x = numcell(1:5,:);
x = x(~isnan(x));
y = meannet(1:5,:);
y = y(~isnan(y));
scatter(x,y,'g','filled')
title('Base')
xlabel('# identified cells')
ylabel('Median of net size')
ax = gca;
% ax.XLim = [0 500];
% ax.XTick = 0:100:500;
ax.YLim = [5 15];
ax.FontSize = ftsz-10;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.LineWidth = lnwd;



subplot(1,3,2)
x = numcell(6:8,:);
x = x(~isnan(x));
y = meannet(6:8,:);
y = y(~isnan(y));
scatter(x,y,'g','filled')
title('Cues')
xlabel('# identified cells')
ylabel('Median of net size')
ax = gca;
% ax.XLim = [0 500];
ax.YLim = [5 15];
ax.FontSize = ftsz-10;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.LineWidth = lnwd;

subplot(1,3,3)
x = numcell(9:11,:);
x = x(~isnan(x));
y = meannet(9:11,:);
y = y(~isnan(y));
scatter(x,y,'g','filled')
title('Airpuff')
xlabel('# identified cells')
ylabel('Median of net size')
ax = gca;
% ax.XLim = [0 500];
ax.YLim = [5 15];
ax.FontSize = ftsz-10;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.LineWidth = lnwd;

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\Network median')
%printpdf('Network median')
%%  number of network events box plot
figure('name','Network Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Fire.mouseID;
yin = Baseline.Network.netoccur(:,2:3);
yin(cell2mat(mouseID(:,2))==1,:) = NaN;
yin(cell2mat(mouseID(:,2))==2,:) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end
% y1 = Cues.Fire.meanfire(cell2mat(Cues.Fire.mouseID(:,2))==trial,2:3);
mouseID = Cues.Fire.mouseID;
yin =  Cues.Network.netoccur(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Fire.mouseID;
yin =  Airpuff.Network.netoccur(:,2:3);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y y1 y2];

h = boxplot(ystat,...
    'BoxStyle','outline',...'Labels',{'Running' 'Resting' 'Running' 'Resting' 'Running' 'Resting'},... 
    'Colors','g');

for j=1:5
    set(h(j,:),'LineWidth',lnwd-1);
 end
set(h(6,:),'LineWidth',lnwd+1);

box off
ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.YLim = [0 250];
ax.LineWidth = lnwd;
% ax.XLim = [0 9];
ax.XTick = [1.5 3.5 5.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
% ax.XTickLabel = {};
ylabel('Events session^-^1','FontSize',ftsz-2)
       

exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\net tot')

%%  frequency of network events running and resting boxplot
figure('name','Network Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   500   420])
hold on
x = [1 2 4 5 7 8];

mouseID = Baseline.Fire.mouseID;
yin = Baseline.Network.netfreq(:,2:3);
yin(cell2mat(mouseID(:,2))==1,:) = NaN;
yin(cell2mat(mouseID(:,2))==2,:) = NaN;
y = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,4))==i,:),1);
end
% y1 = Cues.Fire.meanfire(cell2mat(Cues.Fire.mouseID(:,2))==trial,2:3);
mouseID = Cues.Fire.mouseID;
yin =  Cues.Network.netfreq(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,4))==i,:),1);
end

mouseID = Airpuff.Fire.mouseID;
yin =  Airpuff.Network.netfreq(:,2:3);
y2 = zeros(max(cell2mat(mouseID(:,4))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,4)))
    y2(i,:) = mean(yin(cell2mat(mouseID(:,4))==i & cell2mat(mouseID(:,2))<4,:),1);
end

ystat = [y y1 y2];

h = boxplot(ystat,...
    'BoxStyle','outline',...'Labels',{'Running' 'Resting' 'Running' 'Resting' 'Running' 'Resting'},... 
    'Colors','g');

for j=1:5
    set(h(j,:),'LineWidth',lnwd-1);
 end
set(h(6,:),'LineWidth',lnwd+1);

box off
ax = gca;
ax.FontSize = ftsz-2;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.YLim = [0 15];
ax.LineWidth = lnwd;
% ax.XLim = [0 9];
ax.XTick = [1.5 3.5 5.5];
ax.XTickLabel = {'Baseline' 'Cues' 'AP'};
% ax.XTickLabel = {};
ylabel('Events min^-^1','FontSize',ftsz-2)
       

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\net freq')
%% probality of participations
figure('name','Number of Network events per cell','position',[550   150   950   820])

histogram(Network.netsum,'BinWidth',1,...
    'Normalization','probability',...%     
    'BinLimits',[1, 30],...
    'FaceColor','g',...
    'FaceAlpha',1,...
    'Normalization','probability')
set(gca,'Xlim',[1 30],...  
            'fontsize',35,...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Color',[0 0 0],...
            'LineWidth',2);
% set(gca,'XTick',0:5:maxcount,...
%             'XTickLabel',xtil);
set(gcf,'color',[0 0 0])
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\netprob')


%% Cells per network event
figure('name','Number of cells per Network events','position',[550   150   950   820])

histogram(Network.netnum,'BinWidth',1,...
    'Normalization','probability',...%     
    'BinLimits',[1, 30],...
    'FaceColor','g',...
    'FaceAlpha',1,...
    'Normalization','probability')
set(gca,'Xlim',[1 30],...  
            'fontsize',35,...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Color',[0 0 0],...
            'LineWidth',2);
% set(gca,'XTick',0:5:maxcount,...
%             'XTickLabel',xtil);
set(gcf,'color',[0 0 0])

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\netnum')

%% Box Plot Bulk Baseline

figure('name','Basline Bulksignal',...
    'position',[550   150   650   820],...
    'color',[0 0 0])
hold on

for i = 1:size(Bulk.bulkbase,1)
    plot([1 2],[Bulk.bulkbase(i,3) Bulk.bulkbase(i,2)],...
            'color',[.5 .5 .5],...
            'LineWidth',2)
end

h = boxplot([Bulk.bulkbase(:,3) Bulk.bulkbase(:,2)],...
    'BoxStyle','outline',...
    'Colors','r',...
    'Labels',{'Running' 'Resting'});

 for j=1:5
    set(h(j,:),'LineWidth',3);
 end
set(h(6,:),'LineWidth',6);
set(gca,'fontsize',35,...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Color',[0 0 0],...
            'LineWidth',2)
        

hold off 
[h,p] = ttest(Bulk.bulkbase(:,2),Bulk.bulkbase(:,3));
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\bulkbase')

%% Box plot Corr coefficient
figure('position',[680   400   500   570])
h = boxplot(Bulk.bulkspeed(Bulk.bulkspeed(:,2)<.05,1),...
    'BoxStyle','outline',...
    'Colors','r',...
    'Labels',{''});
 for j=1:5
    set(h(j,:),'LineWidth',3);
 end
 set(h(6,:),'LineWidth',6);
set(gca,'fontsize',35,...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Color',[0 0 0],...
            'LineWidth',2)
        
set(gcf,'color',[0 0 0])

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\bulkruncorr')
%% Delay of bulk signal to running

figure('position',[680   400   500   570])
h = boxplot(Bulk.bulkspeed(Bulk.bulkspeed(:,2)<.05,3),...
    'BoxStyle','outline',...
    'Colors','r',...
    'Labels',{''});
 for j=1:5
    set(h(j,:),'LineWidth',3);
 end
 set(h(6,:),'LineWidth',6);
set(gca,'fontsize',35,...
            'Ylim',[-1 1],...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Color',[0 0 0],...
            'LineWidth',2)
        
set(gcf,'color',[0 0 0])

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\bulkrundelay')

%% bulk network correlation
figure('position',[680   400   400   570])
h = boxplot(Baseline.Bulk.bulknet(Baseline.Bulk.bulknet(:,2)<.05,1),...
    'BoxStyle','outline',...
    'Colors',[1 1 .1],...
    'Labels',{''});
 for j=1:5
    set(h(j,:),'LineWidth',3);
 end
set(h(6,:),'LineWidth',6);
set(gca,'fontsize',35,...
            'YLim',[0 0.25],...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Color',[0 0 0],...
            'LineWidth',2)
        
set(gcf,'color',[0 0 0])

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\bulknetcorr')
%% Delay of summed somatic signal to input bulk signal
figure('position',[680   400   400   570])
h = boxplot(Baseline.Bulk.bulknet(Baseline.Bulk.bulknet(:,2)<.05,3),...
    'BoxStyle','outline',...
    'Colors',[1 1 .1],...
    'Labels',{''});
 for j=1:5
    set(h(j,:),'LineWidth',3);
 end
 set(h(6,:),'LineWidth',6);
set(gca,'fontsize',35,...            
            'YLim',[-.1 .25],...
            'YColor',[1 1 1],...
            'XColor',[1 1 1],...
            'Color',[0 0 0],...
            'LineWidth',2)
        
set(gcf,'color',[0 0 0])
% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\bulknetdelay')





