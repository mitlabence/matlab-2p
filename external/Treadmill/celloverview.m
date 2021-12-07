%addpath('C:\Users\Admin\Documents\MATLAB\export fig')
mouse = [1 2 3 4 5 6 7 8 9];
% mouse = [1 4 5 7 8 ];
exp = [3 4 5 6 7 8];
mouseID = {'M103', 'M155','M158','M194','M195','M224','M226','M227','M234'}; %'M229',
expID = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};

remap = nan(3,length(exp),length(mouse));
corrpossum = nan(3,length(exp),length(mouse));
corrposmean = nan(3,length(exp),length(mouse));
corrplacepossum = nan(3,length(exp),length(mouse));
corrplaceposmean = nan(3,length(exp),length(mouse));
corrstimpossum = nan(3,length(exp),length(mouse));
corrstimposmean = nan(3,length(exp),length(mouse));
daycorrplacepossum = nan(3,length(exp),length(mouse));
daycorrplaceposmean = nan(3,length(exp),length(mouse));
partnerhist = [];
ClustProp = cell(0,0);

for i = 2:length(mouse)
    
%     figure('color',[1 1 1],...
%         'position',[300 100 1.5*[ 594 420]],...
%         'visible','on',...
%         'PaperUnits','centimeters',...
%         'Units','centimeters')
%     
%     cellraster(CAIMcorr,mouse(i),exp);
%     title(mouseID{mouse(i)});

%     cellindi(CAIMcorr,mouse(i),exp);  
    
   remapOUT = cellremap(CAIMcorr,mouse(i),exp);
   remap(:,:,i) = remapOUT;
   
   [cellcorr,cellDay] = cellcorrdays(CAIMcorr,mouse(i),exp);
%    cellcorr = cellDay;
   corrpossum(:,1:size(cellcorr.corrpossum,2),i) = cellcorr.corrpossum;
   corrposmean(:,1:size(cellcorr.corrpossum,2),i) = cellcorr.corrposmean;
   corrplacepossum(:,1:size(cellcorr.corrpossum,2),i) = cellcorr.corrplacepossum;
   corrplaceposmean(:,1:size(cellcorr.corrpossum,2),i) = cellcorr.corrplaceposmean;
   corrstimpossum(:,1:size(cellcorr.corrpossum,2),i) = cellcorr.corrstimpossum;
   corrstimposmean(:,1:size(cellcorr.corrpossum,2),i) = cellcorr.corrstimposmean;
   daycorrplacepossum(:,1:size(cellcorr.corrpossum,2),i) = cellDay.corrplacepossum;
   daycorrplaceposmean(:,1:size(cellcorr.corrpossum,2),i) = cellDay.corrplaceposmean;
   
%    cluster idendity

   [ClustProptemp,partnerhisttemp] = cellclusterdays(CAIMcorr,mouse(i),exp);
%    ClustProp = [ClustProp ClustProptemp];
%    partnerhist = [partnerhist; partnerhisttemp];
end

%%
% files = dir('*.pdf');
% append_pdfs(['PlaceCells.pdf'],files.name)
% delete(files.name)


%% Place field remapping overview
remaptemp = nansum(remap,3);
figure('color',[1 1 1],...
    'position',[500 50 1.5*[420 594]],...
    'renderer','painters',...
    'visible','on')
subplot(3,1,1)
title('Fields of PCs')
hold on
plot(remaptemp(1,:),'color',[.5 .5 .5]);
plot(remaptemp(2,:),'r');
plot(remaptemp(3,:),'g');
ax = gca;
ax.XTick  = 1:size(remaptemp,2);
ax.XTickLabel = expID(exp);
% boxplot(remaptemp,...
%     'BoxStyle','outline',...
%     'Colors','b',...
%     'Labels',{'First field' 'Same field' 'New field'});
ylabel('Mean number of fields');
legend('new PCs','same field','changed field')
grid on

subplot(6,1,3)
title('Followed place cells to place cells')
hold on
y = permute(corrplaceposmean(2,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ylabel('mean r')
ax = gca;
ax.YLim = [0 1];
ax.XTickLabel = {};

subplot(6,1,4)
% title('Followed place cells to place cells')
hold on
y = permute(corrplacepossum(2,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 6];
ax.XTickLabel = {};
ylabel('mean #partner')

subplot(6,1,5)
title('Daily place cells to place cells')
hold on
y = permute(daycorrplaceposmean(2,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ylabel('mean r')
ax = gca;
ax.YLim = [0 1];
ax.XTickLabel = {};

subplot(6,1,6)
% title('Daily place cells to place cells')
hold on
y = permute(daycorrplacepossum(2,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));

ax = gca;
ax.YLim = [0 6];
ax.XTickLabel = {};
ylabel('mean #partner')

% printpdf('Remap')

%% Correlations summary figure
figure('color',[1 1 1],...
    'position',[500 50 1.5*[420 594]],...
    'renderer','painters',...
    'visible','on')

subplot(6,3,1)
title('Stim cells to stim cells')
hold on
y = permute(corrstimposmean(3,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};
ylabel('Mean R');

subplot(6,3,2)
title('Place cells to stim cells')
hold on
y = permute(corrstimposmean(2,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};

subplot(6,3,3)
title('All cells to stim cells')
hold on
y = permute(corrstimposmean(1,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};

subplot(6,3,4)
hold on
y = permute(corrstimpossum(3,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 10];
ylabel('Mean #');

subplot(6,3,5)
hold on
y = permute(corrstimpossum(2,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 10];
% ylabel('Mean #');

subplot(6,3,6)
hold on
y = permute(corrstimpossum(1,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 10];
% ylabel('Mean #');


subplot(6,3,7)
title('Stim cells to place cells')
hold on
y = permute(corrplaceposmean(3,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};
ylabel('Mean R');

subplot(6,3,8)
title('Place cells to place cells')
hold on
y = permute(corrplaceposmean(2,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 .8];
ax.XTickLabel = {};

subplot(6,3,9)
title('All cells to place cells')
hold on
y = permute(corrplaceposmean(1,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};

subplot(6,3,10)
hold on
y = permute(corrplacepossum(3,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 5];
ylabel('Mean #');

subplot(6,3,11)
hold on
y = permute(corrplacepossum(2,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 5];
% ylabel('Mean #');

subplot(6,3,12)
hold on
y = permute(corrplacepossum(1,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 5];
% ylabel('Mean #');

subplot(6,3,13)
title('Stim cells to all cells')
hold on
y = permute(corrposmean(3,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...'Labels',{'Air 1' 'Air 2' 'Air 3'}
    'Colors','b');

ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};
ylabel('Mean R');

subplot(6,3,14)
title('Place cells to all cells')
hold on
y = permute(corrposmean(2,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...'Labels',{'Air 1' 'Air 2' 'Air 3'}
    'Colors','b');

ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};

subplot(6,3,15)
title('All cells to all cells')
hold on
y = permute(corrposmean(1,:,:),[2 3 1]);
% yy = nanmean(permute(corrstimpossum(:,:,:),[2 1 3]),3);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',... 'Labels',{'Air 1' 'Air 2' 'Air 3'}
    'Colors','b');  

ax = gca;
ax.YLim = [0 .5];
ax.XTickLabel = {};

subplot(6,3,16)
hold on
y = permute(corrpossum(3,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 60];
ylabel('Mean #');

subplot(6,3,17)
hold on
y = permute(corrpossum(2,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 60];
% ylabel('Mean #');

subplot(6,3,18)
hold on
y = permute(corrpossum(1,:,:),[2 3 1]);
plot(y,'color',[.5 .5 .5])
boxplot(y',...
    'BoxStyle','outline',...
    'Colors','b',...
    'Labels',expID(exp));
ax = gca;
ax.YLim = [0 60];
% ylabel('Mean #');

% % printpdf('Correlated Cells in AP experiements')
