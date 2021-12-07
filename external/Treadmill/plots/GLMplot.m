
experiment = {'baseline','ICA'};
mouseID = { '177' '184' '235' '255' '339' '239' '342' '349'};
load('/media/2Photon/Nicola/Analisi2020/BigFatSummary.mat')
load('/media/2Photon/Nicola/Analisi2020/BigFatPCA.mat')
%% Performance of linear model learning

num = [5 6]; % decide sessions
MoS = 1;     % Mean(1) or Std(2)
int = 4;     % wihich interval to to test 1:3 learning all/run/rest 4:6 testing all/run/rest
glmp = NaN(5,length(num),2,size(CAIM,2));
mdlperf = NaN(5,length(num),2,size(CAIM,2));
numrounds = NaN(length(num),size(CAIM,2));
numcells = NaN(length(num),size(CAIM,2));
compperf = zeros(5,20);
for i = 1:length(num)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).network)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
%     figure('color',[1 1 1],...
%         'renderer','painters',...
%         'visible','on',...
%         'Units','centimeters',...
%         'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
%         'PaperUnits','centimeters',...
%         'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])
    
    fullCAIM(emptyCAIM) = [];    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        numrounds(i,k) = CAIM(num(i),k).behave.numrounds; 
        numcells(i,k) = size(CAIM(num(i),k).A,2);
        GLMout = PCA(num(i),k).GLMout;
        if ~isempty(GLMout)
            
            y = GLMout.glmp;
            % 1.var 2.ExpID 3.mean or std 4.mouseID
            glmp(:,i,:,k) = y(1:5,int,:); 
            glmp([1:2 5],i,:,k) = y([1:2 5],int+1,:); 
            
            y  = GLMout.mdlperf;
            % 1.var 2.ExpID 3.mean or std 4.mouseID
            mdlperf(:,i,:,k) = y(1:5,int,:);
            % dependence to number of comps
%             compnum = linspace(log(2),log(size(CAIM(num(i),k).cclust,1)),10);
%             compnum = round(exp(compnum));
            compnum = GLMout.numcomp;
            y = GLMout.compperf;
            y = y(:,int,2,:);
            y(4,:) = GLMout.compperf(4,int-1,1,:);
%             y = y(1:100);
%             y = y./y(:,1,1,end);
%             y = y./min(y,[],4);
%             y = y./max(y,[],4);
            compperf(:,:) = y;
            
%             for kk = 1:5
%                 subplot(5,1,kk)
%                 plot(compnum,compperf(kk,:));
%                 hold on
%                 ax = gca;
% %                 ax.YLim = [0 1];
%             end           
        end
    end
end


%% p-value plot
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

subplot(3,3,1)
y = permute(glmp(1,:,2,:),[4 2 1 3]);
boxplot(y)
title(['p of sin, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
ylim([-.05 1.05]);
ax = gca;
ax.XTickLabel = experiment(1:2);

subplot(3,3,2)
y = permute(glmp(2,:,2,:),[4 2 1 3]);
boxplot(y)
title(['p of cos, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
ylim([-.05 1.05]);
ax = gca;
ax.XTickLabel = experiment(1:2);

subplot(3,3,3)
y = permute(glmp(5,:,2,:),[4 2 1 3]);
boxplot(y)
title(['p of theta, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
ylim([-.05 1.05]);
ax = gca;
ax.XTickLabel = experiment(1:2);

subplot(3,2,3)
y = permute(glmp(3,:,2,:),[4 2 1 3]);
boxplot(y)
title(['p of speed, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
ylim([-.05 1.05]);
ax = gca;
ax.XTickLabel = experiment(1:2);

subplot(3,2,4)
y = permute(glmp(4,:,1,:),[4 2 1 3]);
boxplot(y)
title(['p of classifier, ' num2str(round(100*sum(y<=.05)/sum(~isnan(y(:,1))))) ' % sig'])
ylim([-.05 1.05]);
ax = gca;
ax.XTickLabel = experiment(1:2);

subplot(3,2,5)
y = permute(glmp(5,:,2,:),[2 4 1 3]);
scatter(numrounds(:),y(:))
ax = gca;
xlabel('rounds/sessio')

subplot(3,2,6)
y = permute(glmp(5,:,2,:),[2 4 1 3]);
scatter(numcells(:),y(:))
ax = gca;
xlabel('cells/session')

% print(gcf, '-dpdf','GLM p-value')
%% mean of glm
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

a = zeros(1,size(CAIM,2));
b = zeros(1,size(CAIM,2));
MoS = 1;
subplot(5,1,1)
a(:,:) = mdlperf(1,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(1,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('mean of sin')

subplot(5,1,2)
a(:,:) = mdlperf(2,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(2,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('mean of cos')
subplot(5,1,3)
a(:,:) = mdlperf(3,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(3,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('mean of speed')
subplot(5,1,4)
a(:,:) = mdlperf(4,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(4,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('mean of classifier')
subplot(5,1,5)
a(:,:) = mdlperf(5,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(5,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('mean of position')

% print(gcf, '-dpdf','GLM performance mean')
%% std of glm
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

a = zeros(1,size(CAIM,2));
b = zeros(1,size(CAIM,2));
MoS = 2;
subplot(5,1,1)
a(:,:) = mdlperf(1,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(1,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('std of sin')
subplot(5,1,2)
a(:,:) = mdlperf(2,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(2,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('std of cos')
subplot(5,1,3)
a(:,:) = mdlperf(3,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(3,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('std of speed')
subplot(5,1,4)
a(:,:) = mdlperf(4,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(4,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('std of classifier')
subplot(5,1,5)
a(:,:) = mdlperf(5,1,MoS,:);
boxplot(a);hold on
b(:,:) = mdlperf(5,2,MoS,:);
scatter(1:size(CAIM,2),b,'g');hold off
title('std of position')

% print(gcf, '-dpdf','GLM performance std')
%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])
a = zeros(2,size(CAIM,2));
mdlperf1 = cat(2,nanmean(mdlperf(:,1,:,:),2),nanmean(mdlperf(:,2,:,:),2));

subplot(5,2,1)
a(:,:) = mdlperf1(1,:,1,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('mean of sin')
subplot(5,2,2)
a(:,:) = mdlperf1(1,:,2,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('std of sin')
subplot(5,2,3)
a(:,:) = mdlperf1(2,:,1,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('mean of cos')
subplot(5,2,4)
a(:,:) = mdlperf1(2,:,2,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('std of cos')
subplot(5,2,5)
a(:,:) = mdlperf1(3,:,1,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('mean of speed')
subplot(5,2,6)
a(:,:) = mdlperf1(3,:,2,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('std of speed')
subplot(5,2,7)
a(:,:) = mdlperf1(4,:,1,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('mean of class')
subplot(5,2,8)
a(:,:) = mdlperf1(4,:,2,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('std of class')
subplot(5,2,9)
a(:,:) = mdlperf1(5,:,1,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('mean of pos')
subplot(5,2,10)
a(:,:) = mdlperf1(5,:,2,:);
boxplot(a')
hold on
plot(a,'color',[.5 .5 .5])
title('std of pos')

% print(gcf, '-dpdf','GLM base - Cue')



%% individual component number and shuffle
% mouse = {'CA1 1' 'CA1 2'}; 
int = 5;
% 1.var 2.ExpID 3.mean/std 4.mouseID

Readout = 6:10;

num = [5 6];

Position = nan(15,length(num),size(CAIM,2));
Speed = nan(15,length(num),size(CAIM,2));
Classi = nan(15,length(num),size(CAIM,2));

for i = 1:length(num)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).network) || isempty(CAIM(num(i),fullCAIM(j)).plcfield)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end

    fullCAIM(emptyCAIM) = [];
    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);

        x = PCA(num(i),k).GLMout.numcomp;
        if length(x)<15;x(end+1:15) = nan;end
        y = PCA(num(i),k).GLMout.compperf;
        MoS = 2;
        y = y(Readout,int,MoS,:);
             y(4,1,1,:) = PCA(num(i),k).GLMout.compperf(4,4,1,:);
        y = permute(y,[1 4 2 3]);
%         yy = PCA(num(i),k).GLMout.compshuffleMean;
%         yy = yy(Readout,int,MoS,:,:);
%             yy(4,1,1,:,:) = PCA(num(i),k).GLMout.compshuffleMean(4,4,1,:);
%         yyerr = std(yy,0,5);
%         yyerr = permute(yyerr,[1 4 2 3]);
% %         yy = mean(yy,5);
%         yy = permute(yy,[1 4 2 3]);
        %%
        if size(y,2)<15;y(:,end+1:15)= nan;end
        Speed(:,i,k)    = y(3,:);%(y(3,:)/min(y(3,:)));
        Classi(:,i,k)   = y(4,:);%(y(4,:)/min(y(4,:)));
        Position(:,i,k) = y(5,:)./pi*150;%(y(5,:);/min(y(5,:))
%         Position(:,i,k) = mat2gray(Position(:,i,k));
        %%
%         figure('color',[1 1 1],...
%                 'renderer','painters',...
%                 'visible','on',...
%                 'Units','centimeters',...
%                 'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
%                 'PaperUnits','centimeters',...
%                 'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
%             
%         for ii = 1:length(Readout)            
%             subplot(5,1,ii)
%             hold on
%             plot(x,y(ii,:))
% %             errorbar(x,yy(ii,:),yyerr(ii,:))
%         end
%         
%         subplot(5,1,1)
%         title([mouseID(k) ' - sin position'])
%         ylabel('std(predict(sin)-real(sin))')
%         xlabel('number of components')
%         subplot(5,1,2)
%         title('cos position')
%         ylabel('std(predict(cos)-real(cos))')
%         xlabel('number of components')
%         subplot(5,1,3)
%         title('Speed')
%         ylabel('std(predict(speed)-real(speed))')
%         xlabel('number of components')
%         subplot(5,1,4)
%         title('Predictor')
%         ylabel('% false predictions')
%         xlabel('number of components')
%         subplot(5,1,5)
%         title('position')
%         ylabel('std(predict(ang)-real(ang))')
%         xlabel('number of components')
%         print(gcf, '-dpdf',['GLM performance ' mouseID{k}])
    end
end
%%
 figure('color',[1 1 1],...
                'renderer','painters',...
                'visible','on',...
                'Units','centimeters',...
                'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
                'PaperUnits','centimeters',...
                'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
% Position = Classi;
subplot(2,1,1)
plot(x,nanmean(Position(:,1,1:2),[3 2]),'color',[.8 .1 .5])
hold on
plot(x,nanmean(Position(:,1,3:6),[3 2]),'color',[.3 .78 .2])
plot(x,nanmean(Position(:,2,3:6),[3 2]),'color',[.3 .08 .67])
hold off
legend({'WT - Base' 'KA - Base' 'KA - ICA'})

subplot(2,1,2)
for i = 1:size(Position,2)-1
    plot(x,permute(Position(:,i,1:2),[1 3 2]),'color',[.8 .1 .5])
    hold on
    plot(x,permute(Position(:,i,3:6),[1 3 2]),'color',[.3 .78 .2])
end

plot(x,permute(Position(:,2,3:6),[1 3 2]),'color',[.3 .08 .67])

hold off

 print(gcf, '-dpdf','Cell num')
%% Comparisson Baseline Control vs. KA

% mouse = {'M103' 'M155' 'M158' 'M194' 'M195' 'M224' 'M226' 'M227' 'M229' 'M234'}; 
% mouse = {'CA1 1' 'CA1 2'}; 
int = 5;
% 1.var 2.ExpID 3.mean/std 4.mouseID
MoS = 2;
Readout = 6:10;

num = [5];
perfBs = [];
perferrBs = [];
for i = 1:length(num)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).A) || (CAIM(num(i),fullCAIM(j)).behave.numrounds)<5
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end

    fullCAIM(emptyCAIM) = [];    
    for j = 1:2%3:length(fullCAIM)    
        k = fullCAIM(j);

        y = PCA(num(i),k).GLMout.mdlperf;
        y = y(Readout,int,MoS);
            y(4,1,1) = PCA(num(i),k).GLMout.mdlperf(Readout(4),4,1);
        yy = PCA(num(i),k).GLMout.glmshuffleMean;
        yy = yy(Readout,int,MoS);
            yy(4,1,1) = PCA(num(i),k).GLMout.glmshuffleMean(Readout(4),4,1);
        yyerr = PCA(num(i),k).GLMout.glmshuffleStd;
        yyerr = yyerr(Readout,int,MoS);
            yyerr(4,1,1) = PCA(num(i),k).GLMout.glmshuffleStd(Readout(4),4,1);
        perfBs = [perfBs y];
        perferrBs = [perferrBs yy];
        
    end
end

num = [5];
perfCue = [];
perferrCue = [];
for i = 1:length(num)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).A) || (CAIM(num(i),fullCAIM(j)).behave.numrounds)<5
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end

    fullCAIM(emptyCAIM) = [];    
    for j = 3:length(fullCAIM)    
        k = fullCAIM(j);

        y = PCA(num(i),k).GLMout.mdlperf;
        y = y(Readout,int,MoS);
            y(4,1,1) = PCA(num(i),k).GLMout.mdlperf(Readout(4),4,1);
        yy = PCA(num(i),k).GLMout.glmshuffleMean;
        yy = yy(Readout,int,MoS);
            yy(4,1,1) = PCA(num(i),k).GLMout.glmshuffleMean(Readout(4),4,1);
        yyerr = PCA(num(i),k).GLMout.glmshuffleStd;
        yyerr = yyerr(Readout,int,MoS);
            yyerr(4,1,1) = PCA(num(i),k).GLMout.glmshuffleStd(Readout(4),4,1);
        perfCue = [perfCue y];
        perferrCue = [perferrCue yy];
        
    end
end

%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 10 [ 4*8.9 sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [4*8.9 sqrt(2)*8.9]);

% y = [mean(perfBs,2) mean(perferrBs,2) mean(perfCue,2) mean(perferrCue,2)];
% yerr = [std(perfBs,[],2) std(perferrBs,[],2) std(perfCue,[],2) std(perferrCue,[],2)]/sqrt(size(perfCue,2));

% subplot(1,3,1)
% i = 3;
% y = [perfBs(i,:) perferrBs(i,:) perfCue(i,:) perferrCue(i,:)]';
% group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfBs,2)) 3*ones(1,size(perfCue,2)) 4*ones(1,size(perfCue,2))] ; 
% boxplot(y,group)
% ax = gca;
% ax.XTickLabel = {'Ctr' 'shuffle' 'KA' 'shuffle'};
% ylabel('mean error (cm/s)')
% title('Error Speed')
% 
% subplot(1,3,2)
% i = 5;
% y = [perfBs(i,:) perferrBs(i,:) perfCue(i,:) perferrCue(i,:)]'/pi*150;
% group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfBs,2)) 3*ones(1,size(perfCue,2)) 4*ones(1,size(perfCue,2))] ; 
% boxplot(y,group)
% ax = gca;
% ax.XTickLabel = {'Ctr' 'shuffle' 'KA' 'shuffle'};
% ylabel('mean error postion (cm)')
% title('Error Position')
% 
% subplot(1,3,3)
% i = 4;
% y = [perfBs(i,:) perferrBs(i,:) perfCue(i,:) perferrCue(i,:)]';
% group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfBs,2)) 3*ones(1,size(perfCue,2)) 4*ones(1,size(perfCue,2))] ; 
% boxplot(y,group)
% ax = gca;
% ax.XTickLabel = {'Ctr' 'shuffle' 'KA' 'shuffle'};
% title('Classifier run/rest')

subplot(1,3,1)
i = 3;
y = [perfBs(i,:)  perfCue(i,:)]';
group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfCue,2))] ; 
boxplot(y,group)
ax = gca;
ax.XTickLabel = {'Ctr' 'KA' };
ylabel('mean error (cm/s)')
title('Error Speed')

subplot(1,3,2)
i = 5;
y = [perfBs(i,:) perfCue(i,:)]'/pi*150;
group = [ones(1,size(perfBs,2))  2*ones(1,size(perfCue,2)) ] ; 
boxplot(y,group)
ax = gca;
ax.XTickLabel = {'Ctr' 'KA'};
ylabel('mean error postion (cm)')
title('Error Position')

subplot(1,3,3)
i = 4;
y = [perfBs(i,:) perfCue(i,:) ]';
group = [ones(1,size(perfBs,2))  2*ones(1,size(perfCue,2))] ; 
boxplot(y,group)
ax = gca;
ax.XTickLabel = {'Ctr' 'KA' };
title('Classifier run/rest')

% print(gcf, '-dpdf',['Ctr vs KA GLM performance'])


%% Comparisson KA Basline vs ICA

% mouse = {'M103' 'M155' 'M158' 'M194' 'M195' 'M224' 'M226' 'M227' 'M229' 'M234'}; 
% mouse = {'CA1 1' 'CA1 2'}; 
int = 5;
% 1.var 2.ExpID 3.mean/std 4.mouseID
MoS = 2;
Readout = 6:10;

num = [5];
perfBs = [];
perferrBs = [];
for i = 1:length(num)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).A) || (CAIM(num(i),fullCAIM(j)).behave.numrounds)<5
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end

    fullCAIM(emptyCAIM) = [];    
    for j = 3:length(fullCAIM)    
        k = fullCAIM(j);

        y = PCA(num(i),k).GLMout.mdlperf;
        y = y(Readout,int,MoS);
            y(4,1,1) = PCA(num(i),k).GLMout.mdlperf(Readout(4),4,1);
        yy = PCA(num(i),k).GLMout.glmshuffleMean;
        yy = yy(Readout,int,MoS);
            yy(4,1,1) = PCA(num(i),k).GLMout.glmshuffleMean(Readout(4),4,1);
        yyerr = PCA(num(i),k).GLMout.glmshuffleStd;
        yyerr = yyerr(Readout,int,MoS);
            yyerr(4,1,1) = PCA(num(i),k).GLMout.glmshuffleStd(Readout(4),4,1);
        perfBs = [perfBs y];
        perferrBs = [perferrBs yy];
        
    end
end

num = [6];
perfCue = [];
perferrCue = [];
for i = 1:length(num)
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).A) || (CAIM(num(i),fullCAIM(j)).behave.numrounds)<5
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end

    fullCAIM(emptyCAIM) = [];    
    for j = 3:length(fullCAIM)    
        k = fullCAIM(j);

        y = PCA(num(i),k).GLMout.mdlperf;
        y = y(Readout,int,MoS);
            y(4,1,1) = PCA(num(i),k).GLMout.mdlperf(Readout(4),4,1);
        yy = PCA(num(i),k).GLMout.glmshuffleMean;
        yy = yy(Readout,int,MoS);
            yy(4,1,1) = PCA(num(i),k).GLMout.glmshuffleMean(Readout(4),4,1);
        yyerr = PCA(num(i),k).GLMout.glmshuffleStd;
        yyerr = yyerr(Readout,int,MoS);
            yyerr(4,1,1) = PCA(num(i),k).GLMout.glmshuffleStd(Readout(4),4,1);
        perfCue = [perfCue y];
        perferrCue = [perferrCue yy];
        
    end
end

%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 10 [ 4*8.9 sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [4*8.9 sqrt(2)*8.9]);

% y = [mean(perfBs,2) mean(perferrBs,2) mean(perfCue,2) mean(perferrCue,2)];
% yerr = [std(perfBs,[],2) std(perferrBs,[],2) std(perfCue,[],2) std(perferrCue,[],2)]/sqrt(size(perfCue,2));

% subplot(1,3,1)
% i = 3;
% y = [perfBs(i,:) perferrBs(i,:) perfCue(i,:) perferrCue(i,:)]';
% group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfBs,2)) 3*ones(1,size(perfCue,2)) 4*ones(1,size(perfCue,2))] ; 
% boxplot(y,group)
% ax = gca;
% ax.XTickLabel = {'Base' 'shuffle' 'ICA' 'shuffle'};
% ylabel('mean error (cm/s)')
% title('Error Speed')
% 
% subplot(1,3,2)
% i = 5;
% y = [perfBs(i,:) perferrBs(i,:) perfCue(i,:) perferrCue(i,:)]'./pi*150;
% group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfBs,2)) 3*ones(1,size(perfCue,2)) 4*ones(1,size(perfCue,2))] ; 
% boxplot(y,group)
% ax = gca;
% ax.XTickLabel = {'Base' 'shuffle' 'ICA' 'shuffle'};
% ylabel('mean error postion (cm)')
% title('Error Position')
% 
% subplot(1,3,3)
% i = 4;
% y = [perfBs(i,:) perferrBs(i,:) perfCue(i,:) perferrCue(i,:)]';
% group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfBs,2)) 3*ones(1,size(perfCue,2)) 4*ones(1,size(perfCue,2))] ; 
% boxplot(y,group)
% ax = gca;
% ax.XTickLabel = {'Base' 'shuffle' 'ICA' 'shuffle'};
% title('Classifier run/rest')

subplot(1,3,1)
i = 3;
y = [perfBs(i,:) perfCue(i,:)]';
group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfCue,2))] ; 
boxplot(y,group)
ax = gca;
ax.XTickLabel = {'Base' 'ICA'};
ylabel('mean error (cm/s)')
title('Error Speed')

subplot(1,3,2)
i = 5;
y = [perfBs(i,:)  perfCue(i,:)]'./pi*150;
group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfCue,2))] ; 
boxplot(y,group)
ax = gca;
ax.XTickLabel = {'Base' 'ICA' };
ylabel('mean error postion (cm)')
title('Error Position')

subplot(1,3,3)
i = 4;
y = [perfBs(i,:) perfCue(i,:)]';
group = [ones(1,size(perfBs,2)) 2*ones(1,size(perfCue,2))] ; 
boxplot(y,group)
ax = gca;
ax.XTickLabel = {'Base' 'ICA' };
title('Classifier run/rest')

print(gcf, '-dpdf',['Base vs ICA GLM performance Overview'])







