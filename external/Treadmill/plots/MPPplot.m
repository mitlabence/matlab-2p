%% MPP plots

% % network frequency and occurence plots

Cues = load('Z:\Martin Pofahl\Cues.mat');
Baseline = load('Z:\Martin Pofahl\Baseline.mat');
Airpuff = load('Z:\Martin Pofahl\Airpuff.mat');
load('Z:\Martin Pofahl\BigFatCluster.mat');

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
%% Bar graph MPP amplitude
subplot(3,2,5)

x = [1 2 4 5 7 8 10 11 13 14 16 17];
mouseID = Baseline.Bulk.mouseID;
yin = Baseline.Bulk.bulkbase(:,2:3);
yin(yin(:,2) == 0,:) =NaN;
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

mouseID = Cues.Bulk.mouseID;
yin =  Cues.Bulk.bulkbase(:,2:3);
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
    b.CData(2*i-1,:) = [1 .2 .2];
    b.CData(2*i,:) = [.5 0 0];
end

ax = gca;
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XTick = 1.5 :3: 16.5;
ax.XTickLabel = {'Base1' 'Base2' 'Base3' 'Cues1' 'Cues2' 'Cues3'};

title('MPP amplitude')
%% Bar graph MPP standart deviation
subplot(3,2,6)

x = [1 2 4 5 7 8 10 11 13 14 16 17];
mouseID = Baseline.Bulk.mouseID;
yin = Baseline.Bulk.bulkstd(:,2:3);
yin(yin(:,2) == 0,:) =NaN;
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

mouseID = Cues.Bulk.mouseID;
yin =  Cues.Bulk.bulkstd(:,2:3);
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
    b.CData(2*i-1,:) = [1 .2 .2];
    b.CData(2*i,:) = [.5 0 0];
end

ax = gca;
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XTick = 1.5 :3: 16.5;
ax.XTickLabel = {'Base1' 'Base2' 'Base3' 'Cues1' 'Cues2' 'Cues3'};

title('MPP standard deviation')

%% MPP speed correlation
figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]],...
        'renderer','painters',...
        'visible','on')

num = 3:8;
bincorr = nan(6,9,3,20);
for i = num%1:size(CAIM,1)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).bulk)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
%     intEvInt = [];

    for j = 1:length(fullCAIM)
        k = fullCAIM(j);    
        bincorr(i,k,:,:) = CAIM(i,k).bulk.bincorr;
        x = permute(bincorr(i,k,1,:),[3 4 1 2]);
        y = permute(bincorr(i,k,2,:),[3 4 1 2]);
        yerr = permute(bincorr(i,k,3,:),[3 4 1 2]);
%         bincorr(i,k,:,:) = bincorr(i,k,:,:)./max(bincorr(i,k,:,:),[],4);

        subplot(length(num)+1,size(CAIM,2)+1 ,(i-num(1))*(size(CAIM,2)+1) +k)        
%         plot(.1:.1:20,cumsum(CAIM(i,k).network.intEvInt))
        errorbar(x,y,yerr)
        if i == 3
            title(mouse(k))
        end
        if j == 1
            ylabel([experiment(i) ', df/f'])
        end
%         if i == 8
%             xlabel('time/s')
%         end
        grid on
        axis tight
    end
    
end

Y= [];
for j = 1:size(CAIM,2)    
    x = permute(nanmean(bincorr(:,j,1,:)),[3 4 1 2]);
    y = permute(nanmean(bincorr(:,j,2,:)),[3 4 1 2]);
    Y=[Y;y];
    yerr = permute(bincorr(i,k,3,:),[3 4 1 2]);
    if find(x,1)
        subplot(length(num)+1,size(CAIM,2)+1 ,length(num)*(1+size(CAIM,2)) +j)
        plot(x,y)    

        if j == 1
            ylabel(['Sum, df/f'])
        end

        xlabel('cm/s')

        grid on
        axis tight
    end
end

for i = 3:8     
    x = permute(nanmean(bincorr(i,:,1,:)),[3 4 1 2]);
    y = permute(nanmean(bincorr(i,:,2,:)),[3 4 1 2]);
    yerr = permute(nanstd(bincorr(i,:,2,:)),[3 4 1 2])/sqrt(3);
    if find(x,1)
        subplot(length(num)+1,size(CAIM,2)+1,(i-num(1))*(1+size(CAIM,2))+size(CAIM,2)+1)
        plot(x,y)
        if i == 3
            title('sum')
        end
        grid on
        axis tight
             
    end
end

subplot(7,10,70)
plot(x,mean(Y))
grid on
axis tight