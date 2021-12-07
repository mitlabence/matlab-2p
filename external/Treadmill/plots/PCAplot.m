%%
experiment = {'baseline','ICA'};
mouseID = { '177' '184' '235' '255' '339' '239' '342' '349'};
load('/media/2Photon/Nicola/Analisi2020/BigFatSummary.mat')
load('/media/2Photon/Nicola/Analisi2020/BigFatPCA.mat')

%% Decay of variance explained
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

alpha = [];
ninety = [];
ninetynorm = [];
num = 1:8;
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
               
        
        PCAout = PCA(i,k).PCAout;
        if ~isempty(PCAout)
            alpha(i-num(1)+1,k) = PCAout.alpha;
            y = PCAout.explained;
%             y = cumsum(y)/sum(y);
            thresh = .9;
            ninety(i-num(1)+1,k) = find(abs(y-thresh)== min(abs(y-thresh)));
            ninetynorm(i-num(1)+1,k) = ninety(i-num(1)+1,k)/length(y);
        
        subplot(6,10,(i-num(1))*10+k)
%         plot(y,'color',[ 0 .25 1])
        loglog(y,'color',[ 1 .25 0])
        legend(['\alpha = ' num2str(alpha(i-num(1)+1,k),2)],'Location','SouthEast')
        if i == 3;title(mouseID(j));end
        
%         subplot(6,10,(i-num(1))*10+10)
% %         plot(y,'color',[ .5 .5 .5]);
%         loglog(y,'color',[ .5 .5 .5])
%         hold on
        end
    end
     
%     ax = gca;    
%     hold off
%     if i == 3;title('pool');end
end


% print(gcf, '-dpdf','variance dbllog plot')

%% Alpha value and threshold plot

alpha(alpha==0) = NaN;
ninety(ninety==0) = NaN;
ninetynorm(ninetynorm==0) = NaN;
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 10 [ 2.5*sqrt(2)*8.9 8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 8.9])

subplot(1,3,1)
plot(alpha,'color',[.5 .5 .5]);
hold on
boxplot(alpha')
title('Alpha')
ax = gca;  
ax.XTickLabel = {'base3' 'base4' 'base5' 'cue1' 'cue2' 'cue3'};
subplot(1,3,2)
plot(ninety,'color',[.5 .5 .5]);
hold on
boxplot(ninety')
title('90 % threshold')
ax = gca;  
ax.XTickLabel = {'base3' 'base4' 'base5' 'cue1' 'cue2' 'cue3'};
subplot(1,3,3)
plot(ninetynorm,'color',[.5 .5 .5]);
hold on
boxplot(ninetynorm')
title('Threshold normalized to comp-num')
ax = gca;  
ax.XTickLabel = {'base3' 'base4' 'base5' 'cue1' 'cue2' 'cue3'};

% print(gcf, '-dpdf','Alpha value')

%% Delay net to run
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

delay = [];
num = 3:8;
bin = 1;
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
               
        
        PCAout = PCA(i,k).PCAout;
        if ~isempty(PCAout)
            x = PCAout.delaybin(1:bin:end-1)+bin/2;
            y = PCAout.delay;
            yy = [];
            for ii = 1:length(y)/bin
                yy(1,ii) = sum(y(bin*(ii-1)+1:bin*ii));
            end
            y = yy;
            delay(i-num(1)+1,k,:) = y;
%             subplot(length(num)+1,size(CAIM,2)+1,(i-num(1))*(size(CAIM,2)+1)+k)

            subplot(length(num)+1,size(CAIM,2)+1,(i-num(1))*(size(CAIM,2)+1)+k)
            bar(x,y)
    %         legend(['\alpha = ' num2str(alpha(i-2,k),2)],'Location','SouthEast')
            if i == 3;title(mouse(k));end
%             xlim([-20 20])  
        end
    end
    y(1,:) = nanmean(delay(i-num(1)+1,:,:),2);
    subplot(length(num)+1,size(CAIM,2)+1,(i-num(1))*(size(CAIM,2)+1)+size(CAIM,2)+1)
    bar(x,y);
%     ax = gca;
%     xlim([-20 20])
%     hold off
    if i == 3;title('mean');end
end

for j = 1:size(CAIM,2)
    y(1,:) = nanmean(delay(:,j,:),1);
    subplot(length(num)+1,size(CAIM,2)+1,(i-num(1)+1)*(size(CAIM,2)+1)+j)
    bar(x,y);
%     xlim([-20 20]) 
    if j == 1;ylabel('mean');end
end

y(1,:) = nanmean(nanmean(delay(:,:,:),1),2);
subplot(length(num)+1,size(CAIM,2)+1,(i-num(1)+1)*(size(CAIM,2)+1)+j+1)
bar(x,y);
% xlim([-20 20])

% print(gcf, '-dpdf','delay net to run')
%%
figure
plot(x,smooth(y,1))
% print(gcf, '-dpdf','delay net to run pool close up')
%% Cumulative activity
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

cumrun = NaN(6,10,17000);
cumnet = NaN(6,10,17000);
num = 3:8;
% x =scn.tsscn/1000/60;
for i = num
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A) %%|| isempty(CAIM(i,fullCAIM(j)).PCAout) || sum(CAIM(i,k).network.netev)<10 || CAIM(i,k).behave.numrounds<5
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
   
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
               
        PCAout = PCA(i,k).PCAout;
        
        if ~isempty(PCAout) && sum(CAIM(i,k).network.netev)>10 && CAIM(i,k).behave.numrounds>5
            y  = PCAout.cumrun;       
            yy = PCAout.cumnet;
            

            subplot(length(num)+1,size(CAIM,2)+1,(i-num(1))*(size(CAIM,2)+1)+k)
            plot(y,'color',[0 .21 1])
            hold on
            plot(yy,'color',[1 .21 0])
            if i == 3;title(mouse(k));end     
            y = resample(y,17000,length(y));
            yy = resample(yy,17000,length(yy));
            cumrun(i-num(1)+1,k,:) = y;
            cumnet(i-num(1)+1,k,:) = yy; 
        end
    end
    
    y(:) = nanmean(cumrun(i-num(1)+1,:,:),2);
    yy(:) = nanmean(cumnet(i-num(1)+1,:,:),2);
    subplot(length(num)+1,size(CAIM,2)+1,(i-num(1))*(size(CAIM,2)+1)+(size(CAIM,2)+1))
    plot(y,'color',[0 .21 1])
    hold on
    plot(yy,'color',[1 .21 0])
%     xlim([-20 20])
%     hold off
    if i == 3;title('mean');end
end

for j = 1:size(CAIM,2)
    y(:) = nanmean(cumrun(:,j,:),1);
    yy(:) = nanmean(cumnet(:,j,:),1);
    subplot(length(num)+1,size(CAIM,2)+1,(i-num(1)+1)*(size(CAIM,2)+1)+j)
    plot(y,'color',[0 .21 1])
    hold on
    plot(yy,'color',[1 .21 0])
    if j == 1;ylabel('mean');end
end


y(:) = nanmean(nanmean(cumrun(:,:,:),1),2);
yy(:) = nanmean(nanmean(cumnet(:,:,:),1),2);
subplot(length(num)+1,size(CAIM,2)+1,(i-num(1)+1)*(size(CAIM,2)+1)+j+1)
plot(y,'color',[0 .21 1])
hold on
plot(yy,'color',[1 .21 0])

% xlim([-20 20])
% print(gcf, '-dpdf','Cumulative net and run act')
%% Cumulative activity of components

num = [5];

place = [];
net = [];
for i = num
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).network) 
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
               
        PCAout = PCA(i,k).PCAout;
        if ~isempty(PCAout) && sum(CAIM(i,k).network.netev)>10 && CAIM(i,k).behave.numrounds>5
            abinrun = PCAout.abinrun;
            abinnet = PCAout.abinnet;
            for ii = 1:size(abinrun,2)
                place = [place; find(abinrun(:,ii))];
                net =  [net; find(abinnet(:,ii))];
            end
            
            
        end
    end
end

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

plot(cumsum(histcounts(place,1:1:18000))/sum(histcounts(place,1:1:18000)));
hold on
plot(cumsum(histcounts(net,1:1:18000))/sum(histcounts(net,1:1:18000)));
hold off
legend({'place' 'net'});
% xlabel('Number of net events')
% ylabel('p-value')
% % title('cos-dist GPFA')
% print(gcf, '-dpdf','Cumulative activity of components Baseline')

%% Shuffle of network times
num = [3:5 8];
p = NaN(length(num),size(PCA,2),2);
close all
kk = 11;
for j = 1:size(PCA,2) 
    
    if ~isempty(PCA(num(1),j).PCAout)
        figure('color',[1 1 1],...
            'renderer','painters',...
            'visible','on',...
            'Units','centimeters',...
            'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
            'PaperUnits','centimeters',...
            'PaperSize', [2*8.9 2*sqrt(2)*8.9]) 
    end
    
    for i = 1:length(num) 
        
        PCAout = PCA(num(i),j).PCAout;
        
        if ~isempty(PCAout) && isfield(PCAout,'sectioned') && ~isempty(PCAout.sectioned) && ~isempty(PCAout.sectioned.Simil) 

            
            Simil = PCAout.sectioned.Simil;
            SimShuf = PCAout.sectioned.SimShuf;
            
            p(i,j,1) = sum(SimShuf(:,1)<Simil(1))/length(SimShuf(:,1));
            p(i,j,2) = sum(SimShuf(:,2)<Simil(2))/length(SimShuf(:,2));
            %% figure
            subplot(length(num),2,2*(i-1)+1)            
            histogram(SimShuf(:,1),'normalization','probability')
            hold on
            plot([Simil(1) Simil(1)], [0 .5])
            if i == 1;title([mouse{j} ', ' num2str(p(i,j,1))]);else;title(['p = ' num2str(p(i,j,1))]);end
            subplot(length(num),2,2*(i-1)+2)            
            histogram(SimShuf(:,2),'normalization','probability')
            hold on
            plot([Simil(2) Simil(2)], [0 .5])
            title(['p = ' num2str(p(i,j,2))])
        end
    end
%     print(gcf, '-dpdf',num2str(kk))
%     kk = kk+1;
end

%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*8.9 2*sqrt(2)*8.9]) 

subplot(2,2,1)
x = ones(size(p,2),1);
boxplot(p(:,:,1)')
hold on
scatter(x,p(1,:,1))
scatter(x+1,p(2,:,1))
scatter(x+2,p(3,:,1))
scatter(x+3,p(4,:,1))
title('cos^2')

subplot(2,2,2)
x = ones(size(p,2),1);
boxplot(p(:,:,2)')
hold on
scatter(x,p(1,:,2))
scatter(x+1,p(2,:,2))
scatter(x+2,p(3,:,2))
scatter(x+3,p(4,:,2))
title('EROS')

% print(gcf, '-dpdf',num2str(kk))
% kk = kk+1;
% 
% files = dir('*.pdf');
% append_pdfs(['NetworkTimeShuffle.pdf'],files.name)
% delete(files.name)
% close all


%% p-value of cos-dist

num = [3:5 8];

x = nan(length(num),size(CAIM,2));
x2 = nan(length(num),size(CAIM,2));
p = nan(length(num),size(CAIM,2));
p2 = nan(length(num),size(CAIM,2));
p3 = nan(length(num),size(CAIM,2));
for i = 1:length(num)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).network) 
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
               
        PCAout = PCA(num(i),k).PCAout.sectioned;
        GPFAout = PCA(num(i),k).GPFAout;
        network = CAIM(num(i),k).network;
        if ~isempty(PCAout) && ~isempty(PCAout.dist) && ~isempty(GPFAout.dist)
            x(i,k) = sum(network.netev);
            x2(i,k) = max(CAIM(num(i),k).behave.numrounds);
%             x(i,k) = size(GPFAout.dist,2);
            p(i,k)  = GPFAout.distp;
            p2(i,k) = PCAout.rankp;
            
            if ~isempty(PCAout.dual.distHist)
                p3(i,k) = PCAout.dual.rankp;                
            end
        end
    end
end
% p(p>0.05) = NaN;
% p(~isnan(p)) = 1;

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2*8.9])
subplot(1,2,1)
% scatter(x(:),p(:))

scatter(x(:),p2(:))
hold on
% semilogy(x(:),p3(:),'.')
hold on
scatter(x(:),p3(:))
xlabel('Number of net events')
ylabel('p-value')
legend({'GPFA' 'PCA' 'Dual PCA'})
title('cos-dist p value against number of networks')

subplot(1,2,2)
% semilogy(x2(:),p3(:),'.')

% scatter(x2(:),p(:))
scatter(x2(:),p2(:))
hold on
scatter(x2(:),p3(:))
xlabel('rounds per session')
ylabel('p-value')
legend({'GPFA' 'PCA' 'Dual PCA'})
title('cos-dist p value against number of rounds')
% print(gcf, '-dpdf','Cos-dist pvalues')

[sum(p <= .05)/sum(~isnan(p)) sum(p2 <= .05)/sum(~isnan(p2)) sum(p3(x2>5) <= .05)/sum(~isnan(p3(x2>5)))]
%% 

num = [3:5 8];

k = 10;
Distcount = NaN(length(num),size(PCA,2),201);
Randdist = NaN(length(num),size(PCA,2),201);
distPerc = NaN(length(num),size(PCA,2),2,2);
Signum = NaN(length(num),size(PCA,2));
DistStr = NaN(200,length(num),size(PCA,2));
MeanDist = NaN(length(num),size(PCA,2),2);

for j = 1:size(PCA,2) 
    
%     if ~isempty(PCA(num(1),j).PCAout)
%     figure('color',[1 1 1],...
%         'renderer','painters',...
%         'visible','on',...
%         'Units','centimeters',...
%         'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
%         'PaperUnits','centimeters',...
%         'PaperSize', [2*8.9 2*sqrt(2)*8.9]) 
%     end
    
    for i = 1:length(num) 
        
        PCAout = PCA(num(i),j).PCAout;
        
        if ~isempty(PCAout) && isfield(PCAout,'sectioned') && ~isempty(PCAout.sectioned) && ~isempty(PCAout.sectioned.dual.dist) 
            PCAout = PCA(num(i),j).PCAout.sectioned.dual;
            dist = PCAout.dist;
            pInd = PCAout.pInd;
            distHist = PCAout.distHist;
            confint = PCAout.confint;
            x = distHist(1,:);
            distcount = distHist(2,:);
            randdist = distHist(3,:);
            Distcount(i,j,:) = distcount;
            Randdist(i,j,:) = randdist;
            distPerc(i,j,:,:) = PCAout.distPerc;
            B = pInd(:,:,2).*dist;
            B(isnan(B)) = 0;
            signum = sum(B(:)~=0);
            Signum(i,j) = signum/length(B(:));
            DistStr(:,i,j) = histcounts(B(B~=0),-1:.01:1);         
            MeanDist(i,j,:) =  PCAout.simil;
            
            %% figure
%             subplot(length(num),3,3*(i-1)+1)
%             if i == 1;title(mouse{j});end
%             hold on
%             plot(x,cumsum(randdist))
%             plot(x,cumsum(distcount))
%             plot([confint(1) confint(1)],[0 1],'color',[.5 .5 .5])
%             plot([-1 1],[.05 .05],'color',[.5 .5 .5])
%             plot([confint(2) confint(2)],[0 1],'color',[.5 .5 .5])
%             plot([-1 1],[.95 .95],'color',[.5 .5 .5])
%             legend({ 'shuf' 'real'},'location','southeast')
%             axis tight
%             xlim([-.5 .5])        
%             
%             subplot(length(num),3,3*(i-1)+2)
%             histogram(B(B~=0),-1:.05:1);
%             
%             b = subplot(length(num),3,3*(i-1)+3);           
% %             B = abs(B);
%             mycolormap = jet(256);
%             mycolormap(1,:) = [1 1 1];
%             imagesc(abs(B));
%             colormap(b,mycolormap);           
%             temp = get(gca,'position');
%             colorbar
%             set(gca,'position',temp)
%             
%             title([num2str(signum) '(' num2str(round(100*signum/length(B(:)))) ' %)' ' sig distances'])

            
        end
    end
%     k = k+1;
%     printpdf([num2str(k)])
end

% y = MeanDist(:,:,2);
% y = Signum;
% group = cell(size(y));
% group(1:3,:) = experiment(5);
% group(4,:) = experiment(8);
% y = [nanmean(y(1:3,:),1); y(4,:)];
% y = y./sum(y,1);
% [h,p] = ttest(y(1,:),y(2,:))
% anova1(y(:),group(:))
% printpdf(['Similarity factor'])

SigBase = Signum(1:3,:);
SigCue = Signum(4,:);
[nanmean(SigBase(:)) nanstd(SigBase(:)) nanstd(SigBase(:))/sqrt(sum(~isnan(SigBase(:))));
nanmean(SigCue) nanstd(SigCue) nanstd(SigCue)/sqrt(sum(~isnan(SigCue(:))))]
%%
x = nansum(DistStr(:,:,:),3);
x = [nansum(x(:,1:3),2) x(:,4)];
x = x./sum(x);
plot(x)

x = (DistStr(:,:,8));
x = [nansum(x(:,1:3),2) x(:,4)];
x = x./sum(x);
% x = permute(x,[2 3 1]);
plot(x)

[nanmean(nanmean(distPerc(:,:,1,1),1),2) nanmean(nanmean(distPerc(:,:,1,2),1),2);
nanmean(nanmean(distPerc(:,:,2,1),1),2) nanmean(nanmean(distPerc(:,:,2,2),1),2)];

% files = dir('*.pdf');
% append_pdfs(['PCA dual individual Distance.pdf'],files.name)
% delete(files.name)
% close all

%%
num = 1:10;
distcount = nanmean(nanmean(Distcount(4,num,:),1),2);
randdist =  nanmean(nanmean(Randdist(4,num,:),1),2);
distcount = permute(distcount,[1 3 2]);
randdist = permute(randdist,[1 3 2]);

figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
        'PaperUnits','centimeters',...
        'PaperSize', [2*8.9 2*sqrt(2)*8.9])

subplot(3,2,1)
histogram(shufdist,-1:.05:1,'normalization','probability')
hold on
histogram(dist,-1:.05:1,'normalization','probability')
legend({ 'shuf' 'real'})

subplot(3,2,2)
hold on
plot(x,cumsum(randdist))
plot(x,cumsum(distcount))
plot([confint(1) confint(1)], [0 1],'--','color',[0 0 0])
plot([confint(1,2) confint(1,2)], [0 1],'--','color',[0 0 0])
plot([confint(2,1) confint(2,1)], [0 1],'--','color',[.5 .5 .5])
plot([confint(2,2) confint(2,2)], [0 1],'--','color',[.5 .5 .5])
legend({ 'shuf' 'real'},'location','southeast')
axis tight

subplot(3,2,3)
plot(x,smooth(distcount-randdist,10))
grid on
% hold on
% plot(cumsum(randdist))

subplot(3,2,4)
hold on
y = cumsum(distcount);
plot(x,cumsum(randdist))
plot(x,y)
plot([confint(1) confint(1)], [0 1],'--','color',[0 0 0])
plot([confint(1,2) confint(1,2)], [0 1],'--','color',[0 0 0])
plot([confint(2,1) confint(2,1)], [0 1],'--','color',[.5 .5 .5])
plot([confint(2,2) confint(2,2)], [0 1],'--','color',[.5 .5 .5])
legend({ 'shuf' 'real'},'location','southeast')
xlim([0.2 0.7])
ylim([y(121) 1])


subplot(3,2,5)
win = find(confint(1,2)==x,1):length(x);
% a = smooth(cumsum(randdist));
% b = smooth(cumsum(distcount));
% plot(x(win),abs(a(win)-b(win)))
hold on
plot(x(win),randdist(win));
plot(x(win),distcount(win));
title([num2str(sum(distcount(win)-randdist(win))) ' % difference'])

subplot(3,2,6)
win = 1:length(x);%find(confint(1,2)==x,1):length(x);
a = (cumsum(randdist));
b = (cumsum(distcount));
plot(x(win),(sign(x).*(a(win)-b(win))))
grid on

%%
figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','on',...
        'Units','centimeters',...
        'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
        'PaperUnits','centimeters',...
        'PaperSize', [2*8.9 2*sqrt(2)*8.9])
    
num = 1:10;
trial = 1:3;
distcount = nanmean(nanmean(Distcount(trial,num,:),1),2);
randdist =  nanmean(nanmean(Randdist(trial,num,:),1),2);
distcount = permute(distcount,[1 3 2]);
randdist = permute(randdist,[1 3 2]);

subplot(3,2,1)
hold on
plot(x,(randdist))
plot(x,(distcount))
legend({ 'shuf' 'real'},'location','southeast')
axis tight

subplot(3,2,2)
hold on
plot(x,cumsum(randdist))
plot(x,cumsum(distcount))
legend({ 'shuf' 'real'},'location','southeast')
axis tight

subplot(3,2,3)
plot(x,smooth(distcount-randdist,10))
grid on
% hold on
% plot(cumsum(randdist))
title('diff distribution')

subplot(3,2,4)
win = 1:length(x);%find(confint(1,2)==x,1):length(x);
a = (cumsum(randdist));
b = (cumsum(distcount));
plot(x(win),(sign(x).*(a(win)-b(win))))
grid on

num = 1:10;
trial = 4;
distcount = nanmean(nanmean(Distcount(trial,num,:),1),2);
randdist =  nanmean(nanmean(Randdist(trial,num,:),1),2);
distcount = permute(distcount,[1 3 2]);
randdist = permute(randdist,[1 3 2]);

 
subplot(3,2,5)
plot(x,smooth(distcount-randdist,10))
grid on
% hold on
% plot(cumsum(randdist))
title('Cue enriched')

subplot(3,2,6)
win = 1:length(x);%find(confint(1,2)==x,1):length(x);
a = (cumsum(randdist));
b = (cumsum(distcount));
plot(x(win),(sign(x).*(a(win)-b(win))))
grid on

% printpdf(['Distance Pool network'])