function PCAout = PCAprint(caim,scn,pathname,filename,method)
addpath('/home/martin/Documents/MATLAB/export fig')
addpath(genpath('/media/2Photon/Matlab code/Martin Source Code/gpfa_v0203'))
kk = 11;
%%
C = caim.C;C = C./caim.Df;     
% C = caim.S;C(caim.S_bin==0)= 0;

clear dat
dat.trialId = 1;
dat.spikes = C;

%%
if strcmp(method,'pca')
    %%
    runIdx = 2;
    xDim = size(C,1)-1;
    binWidth = 1; 
    kernSD = 5;

    % Extract neural trajectories
    result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                        'kernSDList', kernSD,'binWidth',binWidth);
    if isempty(result);PCAout = [];return;end
    [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);%%
    score = result.kern.estParams.L;
    hasSpikesBool = result.hasSpikesBool;
%    score = estParams.Corth;
    a = [];
    for i = 1:size(seqTrain,2)
        a = [a seqTrain(i).xorth];
    end
    clear b

    for i = 1:size(a,1)
        b(i,:) = interp1(1:size(a,2),a(i,:),1:size(a,2)/size(caim.C,2):size(a,2),'spline');
    end

    b(:,end+1:size(caim.C,2)) = 0;
    a = b';
elseif strcmp(method,'ica')
    %%
    disp('I am performing a really nice ICA, bitches')
    runIdx = 2;
    kernSD = 10;
    xDim = 50;%size(C,1)-1;
    yOut = smoother(C, kernSD,1);
    b = rica(yOut',xDim,'IterationLimit',1000);
    % b = rica(a,q,'IterationLimit',1000);
    a = yOut'*b.TransformWeights;
    [~,aa] = sort(var(a),'descend');
    a = a(:,aa);
    score = b.TransformWeights;
    score = score(:,aa);
    hasSpikesBool = true(size(C,1),1);
end

%%
run = find(scn.speed>0);
space = round(scn.distance(run));
% space = round(scn.totdist(run));
space(space<=0) = 1;
speed = scn.speed(run)*100;
speed(speed<0) = 0;
speed(end+1)=speed(end);
speed = smooth(speed,30);
% speed = speed-min(speed);
maxspeed = max(speed);
speed = round(speed/max(speed));

netpos = find(caim.network.netev);

%%

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 2*8.9])

hold on
mycolormap = hsv(max(space));

for i = 1:length(run)-1
   if diff(run(i:i+1))<5
       plot3([a(run(i),2) a(run(i+1),2)],[a(run(i),3) a(run(i+1),3)],[a(run(i),1) a(run(i+1),1)],'color',mycolormap(space(i),:));
    end  
end
view(45,-45)
title([method ' components plotted against space'])
%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;
%%

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 2*8.9])

hold on
win = -20:20;
mycolormap = hsv(length(win));
% mycolormap = hsv(length(win));mycolormap = [mycolormap(floor(end/2)+1:end,:);mycolormap(1:floor(end/2),:)];

% win = -10:5;
% mycolormap = jet(length(win));

netpos1 = netpos(netpos>abs(win(1)) & netpos<size(C,2)-win(end));
netpos1(diff(netpos1)<length(win)) = [];
for i = 1:length(netpos1)
    for j = win
        jj = netpos1(i)+j;
        plot3([a(jj,1) a(jj+1,1)],[a(jj,2) a(jj+1,2)],[a(jj,3) a(jj+1,3)],'color',mycolormap(j-win(1)+1,:));
    end  
end
view(45,-45)
title([method 'components plotted around network'])
%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;
%%
% figure('color',[1 1 1],...
%     'renderer','painters',...
%     'visible','off',...
%     'Units','centimeters',...
%     'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
%     'PaperUnits','centimeters',...
%     'PaperSize', [2*sqrt(2)*8.9 2*8.9])
% 
% hold on
% 
% num = 1:20;
% thresh = 0.02;
% mycolormap = jet(length(num));
% for i = num
%     x = scn.distance(run);
%     y = a(run,i);
%     z = scn.rounds(run);
% %     y = abs(y);
% %     y = rescale(y,0, 1);
% %     gaussEqn = 'a*exp(-((x-b)/(sqrt(2)*c))^2)';
% %     startPoints = [1 x(max(y)==y) 100];
% %     upperBounds = [1 1500 500];
% %     lowerBounds = [0 0 30];
% %     f1 = fit(x,y,gaussEqn,...
% %          'Start', startPoints,... 
% %          'Upper',upperBounds,...  
% %          'Lower',lowerBounds);
% %     plot(1:1500,f1(1:1500));
% %     hold on
%     yy = zeros(1,150);
%     for j = 1:150
%         yy(j) = mean(y(x>(j-1)*10&x<j*10));
%     end
%     x = x((y)>thresh);
%     z = z((y)>thresh);
%     y = y((y)>thresh);
%     scatter3(x,y,z,'filled','markerfacecolor',mycolormap(i-num(1)+1,:))
%     
% %     plot3(x,y,z,'color',mycolormap(i-num(1)+1,:))
% 
% %     scatter(x,y,'filled','markerfacecolor',mycolormap(i-num(1)+1,:))
% %     plot(1:10:1500,yy)
%     
% end
% ylabel('intensity')
% xlabel('distance/mm')
% title('PCA components in different rounds')
% %%
% print(gcf, '-dpdf',num2str(kk))
% kk = kk+1;
%% Components in rounds 

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 2*8.9])

imagesc(scn.totdist(run),1:20,abs(a(run,1:20)'))
ax = gca;
ax.XTick = 1500:1500:max(scn.totdist(run));
ax.XTickLabel = 1:max(scn.totdist(run))/1500;
colormap(jet)
ylabel('PCA ID')
xlabel('round')
title([method ' components in rounds'])
%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;
%% Components trigered to network

netpos = find(caim.network.netev);
span = -15:15;
netnum = 1:length(netpos);
num = 1:size(a,2);
netpos = find(caim.network.netev);
maxcomp = zeros(length(num),length(netnum));
for i = num
    for j = netnum
        if span(end)+netpos(j)>size(C,2)
            int = netpos(j)+span(1):size(C,2);            
        elseif span(1)+netpos(j)<1
            int = 1:netpos(j)+span(end);
        else
            int = span+netpos(j);
        end
        x = span/15;
        y = a(int,i);

        if max(y)>abs(min(y))
            maxcomp(i,j) = max(y);
        elseif isempty(min(y))
            maxcomp(i,j) = 0;
        else
            maxcomp(i,j) = min(y);
        end
    end
end

netpos = find(caim.network.netev);
nettime = caim.network.netpos/1000/60;
if isfield(caim.network,'netraster')
    netraster = caim.network.netraster;
    num = 1:size(a,2);%15;%
    dist = zeros(size(netraster,2),length(num));
    runstart = find(diff(scn.running)==1);
    runstop = find(diff(scn.running)==-1);
    for j = 1:length(num)
        k = num(j);
        for i = 1:size(netraster,2)              
            aa = zeros(size(C,1),1);
            aa(hasSpikesBool) = score(:,k);  
            dist(i,j) = (aa'*netraster(1:end,i))./(sqrt(sum(netraster(:,i)))*sqrt(sum(aa.^2)));
        end
    end
end
%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
% hold on
% for i = num
%     plot(smooth(maxcomp(i,:)),'color',mycolormap(i-num(1)+1,:))
% end
% axis tight
num = 1:30;
if size(a,2)<num(end);num = 1:size(a,2);end

subplot(2,1,1)
imagesc(abs(maxcomp(num,:)))
ylabel('PCA comp ID')
xlabel('Network events')
title(['Intensity of ' method ' in network events'])
colormap(jet)
colorbar
if isfield(caim.network,'netraster')
    subplot(2,1,2)
    imagesc(abs(dist(:,num)'))
    ylabel('PCA comp ID')
    xlabel('Network events')
    title(['Cos-distance between ' method ' components & network events'])
colorbar
end
%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;

%% double logarithmic plot
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
hold on
explained = var(a);

linEqn = 'a+b*(x)';
if length(explained)< 30
    int = 1:round(length(explained));
else
    int = 4:round(length(explained)/3);
end

% x = log(1:length(int));
x = log(int);
y = log(explained(int));
startPoints = [y(1) -1 ];
upperBounds = [y(1)+3 1];
lowerBounds = [y(1)-3 -3 ];
f1 = fit(x',y',linEqn,...
       'Start', startPoints,...
       'Upper',upperBounds,...
       'Lower',lowerBounds);
Alpha = -round(f1.b,2);
int = 1:length(explained);
x = log(1:length(int));
y = log(explained(int));

PCAout.alpha = Alpha;
PCAout.explained = explained;

subplot(2,1,1)
plot(cumsum(explained)/sum(explained))
title(['cumulative sum of variance per component'])
ylabel('var %')
xlabel('# comp')

subplot(2,1,2)
plot(x,y)
hold on
plot(x,f1(x))
hold off
title(['decay of variance per component, 1/x^{-\alpha} fit, \alpha = ' num2str(Alpha)])
ylabel('log(var)/comp')
xlabel('log(# comp)')

% exportfig('fit')
%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;
%%
gaussEqn = 'a*exp(-((x-b)/(sqrt(2)*c))^2)';
bin = max(abs(a(:)))/100;
x = -max(abs(a(:))):bin:max(abs(a(:)));
xx = -(max(abs(a(:))))-bin/2:bin:(max(abs(a(:))))+bin/2;
abin = zeros(size(a));
abinrun = zeros(size(a));
abinnet = zeros(size(a));
adist = NaN([size(a) 2]);
runfirst = false(size(a,2),1);
for i = 1:size(a,2)
    % smoothed actitivity in apsbins for fitting
    y = histcounts(a(:,i),xx,'normalization','probability');
    % maximum as peak strating point
    [a1,b1] = max(y);
    startPoints = [a1 x(b1) 1];
    upperBounds = [1 4*abs(x(b1)) x(end)/4];
    lowerBounds = [0 -4*abs(x(b1)) 0];
    f1 = fit(x',y',gaussEqn,...
        'Start', startPoints,...
        'Upper',upperBounds,...
        'Lower',lowerBounds);
    yy = f1(x);
    abin(a(:,i)>f1.b+3*f1.c,i) = 1;
    abin(a(:,i)<f1.b-3*f1.c,i) = -1; 
    abinrun(abs(abin(:,i))==1 & scn.running==1,i) = abin(abs(abin(:,i))==1 & scn.running==1,i);
    abinnet(abs(abin(:,i))==1 & caim.network.netev==1,i) = abin(abs(abin(:,i))==1 & caim.network.netev==1,i);
    aa = find(abinnet(:,i));
    aa = [1; aa; size(a,1)];
    for j = 2:length(aa)-1
        jj = aa(j);
        if ~isempty(find(abinrun(aa(j-1):jj,i),1,'last'))
            adist(jj,i,1) = (find(abinrun(aa(j-1):jj,i),1,'last')-jj+aa(j-1));
            if j == 2
                runfirst(i) = 1;
            end
        end
        if ~isempty(find(abinrun(jj:aa(j+1),i),1,'first'))
            adist(jj,i,2) = find(abinrun(jj:aa(j+1),i),1,'first');
        end
        if jj>10 && jj< size(a,1)-10
            abinnet(jj-10:jj+10,i)=abinnet(jj,i);
        end
    end
%     figure
% %     histogram(a(scn.running==1,i),xx,'normalization','probability')
% %     hold on
% %     histogram(a(scn.running==0,i),xx,'normalization','probability')
%     plot(x,y)
%     hold on
%     plot(x,yy)
%     plot([f1.b+3*f1.c f1.b+3*f1.c],[0 max(y)])
%     plot([f1.b-3*f1.c f1.b-3*f1.c],[0 max(y)])
%     hold off
end

%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
a1 = zeros(size(a));
a2 = a1;
for i = 1:size(abinrun,2)
    a1(:,i) = smooth(a(:,i).*abs(abinrun(:,i)),1);
    a2(:,i) = smooth(a(:,i).*abs(abinnet(:,i)),1);
end

num = 1:20;
subplot(4,1,1)
imagesc(scn.tsscn/1000/60,num,a1(:,num)',.15*[min(a(:)) max(a(:))])
xlabel('time/min')
ylabel('Component ID')
title('Place Active Components')
colorbar
subplot(4,1,2)
% imagesc(scn.tsscn/1000/60,num,a2(:,num)')
imagesc(scn.tsscn/1000/60,num,a2(:,num)',.05*[min(a(:)) max(a(:))])
xlabel('time/min')
ylabel('Component ID')
title('Network Active Components')

colorbar

clmplg = floor(256*[min(a(:))/(min(a(:))-max(a(:))) 1-min(a(:))/(min(a(:))-max(a(:)))]);


mycolormap = zeros(sum(clmplg),3);
j = clmplg(2);
for i = 1:sum(clmplg)
    if i <clmplg(1)+2
        mycolormap(i,:) = [1 1-(clmplg(1)-(i-1))/clmplg(1) 1-(clmplg(1)-(i-1))/clmplg(1)];  
    else
        j = j-1;
        mycolormap(i,:) = [1-(clmplg(2)-j)/clmplg(2) 1-(clmplg(2)-j)/clmplg(2) 1];
    end
end
        
colormap(mycolormap)
% colormap([1 0 0; 1 1 1; 0 0 1])
% colormap(jet)

a1 = zeros(size(a));
a2 = a1;
for i = 1:size(abinrun,2)
%     a1(:,i) = smooth(a(:,i).*abinrun(:,i),500);
%     a2(:,i) = smooth(abs(a(:,i).*abinnet(:,i)),500);
    a1(:,i) = cumsum(abs(abinrun(:,i)))/sum(abs(abinrun(:,i)));
    a2(:,i) = cumsum(abs(abinnet(:,i)))/sum(abs(abinnet(:,i)));
end

subplot(4,1,3)
plot(scn.tsscn/1000/60,scn.totdist/scn.totdist(end))
hold on
plot(scn.tsscn/1000/60,nanmean(a1,2),'color',[0 .21 1])
plot(scn.tsscn/1000/60,cumsum(caim.network.netev)/sum(caim.network.netev))
plot(scn.tsscn/1000/60,nanmean(a2,2),'color',[1 .21 0])
legend({'Distance run' 'Running components' 'Cumulative network events' 'Network components'},'Location','northwest');
hold off


% subplot(3,3,7)
% histogram(adist(:,num,1),'normalization','probability')
subplot(4,1,4)
histogram(adist(:,:,:)/15,[-200:5:200],'normalization','probability')
title(['\DeltaT net/run, ' 'mean = ' num2str(nanmean(adist(:))/15) ' , ' 'median = ' num2str(nanmedian(adist(:))/15) ' , ' 'ratio(first run/first net) = ' num2str(mean(runfirst))]);
% subplot(3,3,9)
% histogram(adist(:,num,2),'normalization','probability')

%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;

%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*8.9 2*sqrt(2)*8.9]);

for i = 1:5
    subplot(5,5,(i-1)*5+1:(i-1)*5+4)
    runstart = find(diff(scn.running)==1);
    runstop = find(diff(scn.running)==-1);
    
    if isempty(runstart)
        plot(scn.tsscn/1000/60,a(:,i),'color',[1 .21 0])
    else
        if runstop(1)<runstart(1);runstop = runstop(2:end);end
        if runstart(end)>runstop(end);runstart = runstart(1:end-1);end
        runstop = [2;runstop(1:end)];
        hold on

        for j = 1:length(runstart)
            int = runstop(j)-1:runstart(j);
            plot(scn.tsscn(int)/1000/60,a(int,i),'color',[1 .21 0])
            int = runstart(j):runstop(j+1);
            plot(scn.tsscn(int)/1000/60,a(int,i),'color',[0 .21 1])
        end
        int = runstop(j+1):size(a,1);
        plot(scn.tsscn(int)/1000/60,a(int,i),'color',[1 .21 0])
    end
    
    aa = find(abinnet(:,i));
    stpsz = max(a(:,i))/length(aa);
    for j = 1:length(aa)
        jj = aa(j);
        if isnan(adist(jj,i,1)) && ~isnan(adist(jj,i,2))
            plot([scn.tsscn(jj) scn.tsscn(jj + adist(jj,i,2))]/1000/60,[stpsz*j stpsz*j])
        elseif ~isnan(adist(jj,i,1)) && isnan(adist(jj,i,2))
            plot([scn.tsscn(jj+adist(jj,i,1)) scn.tsscn(jj)]/1000/60,[stpsz*j stpsz*j])
        elseif ~isnan(adist(jj,i,1)) && ~isnan(adist(jj,i,2))
            plot([scn.tsscn(jj+adist(jj,i,1)) scn.tsscn(jj + adist(jj,i,2))]/1000/60,[stpsz*j stpsz*j])
        end
        scatter(scn.tsscn(jj)/1000/60,stpsz*j,'filled','markerfacecolor',[0 1 0])
    end
    if i ==1
        title('example traces')
    end
    
%     plot(scn.tsscn/1000/60,-1+scn.totdist/scn.totdist(end))
    hold on
    plot(scn.tsscn/1000/60,-2+2*a1(:,i),'color',[0 .21 1])
%     plot(scn.tsscn/1000/60,-1+cumsum(caim.network.netev)/sum(caim.network.netev))
    plot(scn.tsscn/1000/60,-2+2*a2(:,i),'color',[1 .21 0])
    
    subplot(5,5,(i-1)*5+5)
    histogram(a(scn.running==1 & a(:,i)~=0,i),min(a(:,i)):(max(a(:,i))-min(a(:,i)))/30:max(a(:,i)),'normalization','probability','FaceColor',[0 .21 1])
    hold on
    histogram(a(scn.running==0 & a(:,i)~=0,i),min(a(:,i)):(max(a(:,i))-min(a(:,i)))/30:max(a(:,i)),'normalization','probability','FaceColor',[1 .21 0])
end
%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;

%% plot decoded data (GLM)
if isfield(caim,'GLMout')
    brdln = scn.tsscn(round(4/5*size(caim.C,2)))/1000;
    %% plot traces fitted on raw data
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','off',...
        'Units','centimeters',...
        'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
        'PaperUnits','centimeters',...
        'PaperSize',[2*8.9 2*sqrt(2)*8.9]);
    
    Readout = 1:5;
    
    for ii = 1:length(Readout)           
        subplot(5,1,ii)
        hold on
        y = caim.GLMout.pred(:,Readout(ii));
        plot(scn.tsscn/1000,y)
        yy = caim.GLMout.movement(:,ii);
        plot(scn.tsscn/1000,yy)
        plot([brdln brdln],[min(y) max(y)],'color',[0 .8 .2]);
    end

    subplot(5,1,1)
    title('sin of position, reconstructed vs real, fit from traces')
    ylabel('sin')
    xlabel('time /s')
    subplot(5,1,2)
    title('cos of position, reconstructed vs real')
    ylabel('cos')
    xlabel('time /s')
    subplot(5,1,3)
    title('Speed')
    ylabel('speed, reconstructed vs real')
    xlabel('time /s')
    subplot(5,1,4)
    title('Classifier')
    ylabel('false predictions')
    xlabel('time /s')
    subplot(5,1,5)
    title('position, reconstructed vs real')
    ylabel('ang/\pi')
    xlabel('time /s')

    %%
    print(gcf, '-dpdf',num2str(kk))
    kk = kk+1;

    %% plot error of fit and prediction

    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','off',...
        'Units','centimeters',...
        'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
        'PaperUnits','centimeters',...
        'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
    Readout = 6:10;
    for ii = 1:length(Readout)            
        subplot(5,1,ii)
        hold on
        y = caim.GLMout.test(:,ii);
        plot(scn.tsscn/1000,y)
        plot([brdln brdln],[min(y) max(y)],'color',[.3 .3 .3]);
    end

    subplot(5,1,1)
    title(['error sin position'])
    ylabel('sin(predict)-sin(real)')
    xlabel('number of components')
    subplot(5,1,2)
    title('error cos position')
    ylabel('cos(predict)-cos(real))')
    xlabel('number of components')
    subplot(5,1,3)
    title('error  Speed')
    ylabel('speed(predict)-speed(real)')
    xlabel('number of components')
    subplot(5,1,4)
    title('error Predictor')
    ylabel('false predictions')
    xlabel('number of components')
    subplot(5,1,5)
    title('error position')
    ylabel('ang(predict)-ang(real)')
    xlabel('number of components')

    %%
    print(gcf, '-dpdf',num2str(kk))
    kk = kk+1;
    %% plot traces fitted on PCA traces
%     figure('color',[1 1 1],...
%         'renderer','painters',...
%         'visible','off',...
%         'Units','centimeters',...
%         'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
%         'PaperUnits','centimeters',...
%         'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
%     Readout = 1:5;
%     
%     for ii = 1:length(Readout)           
%         subplot(5,1,ii)
%         hold on
%         y = caim.GLMout.pred(:,Readout(ii));
%         plot(scn.tsscn/1000,y)
%         yy = caim.GLMout.movement(:,ii);
%         plot(scn.tsscn/1000,yy)
%         plot([brdln brdln],[min(y) max(y)],'color',[0 .8 .2]);
%     end
% 
%     subplot(5,1,1)
%     title('sin of position, reconstructed vs real, fit from PCA components')
%     ylabel('sin')
%     xlabel('time /s')
%     subplot(5,1,2)
%     title('cos of position, reconstructed vs real')
%     ylabel('cos')
%     xlabel('time /s')
%     subplot(5,1,3)
%     title('Speed')
%     ylabel('speed, reconstructed vs real')
%     xlabel('time /s')
%     subplot(5,1,4)
%     title('Classifier')
%     ylabel('false predictions')
%     xlabel('time /s')
%     subplot(5,1,5)
%     title('position, reconstructed vs real')
%     ylabel('ang/\pi')
%     xlabel('time /s')
% 
%     %%
%     print(gcf, '-dpdf',num2str(kk))
%     kk = kk+1;
end

%% plot performance against number of components
x = caim.GLMout.numcomp;
y = caim.GLMout.compperf;
Readout = 1:5;
int = 5;
MoS = 2;
y = y(Readout,int,MoS,:);
y(4,1,1,:) = caim.GLMout.compperf(4,4,1,:);
y = permute(y,[1 4 2 3]);
yy = caim.GLMout.compshuffleMean;
yy = yy(Readout,int,MoS,:);
yy(4,1,1,:,:) = caim.GLMout.compshuffleMean(4,4,1,:);
yy = permute(yy,[1 4 2 3]);

yyerr = caim.GLMout.compshuffleStd;
yyerr = yyerr(Readout,int,MoS,:);
yyerr(4,1,1,:,:) = caim.GLMout.compshuffleStd(4,4,1,:);
yyerr = permute(yyerr,[1 4 2 3]);


figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','off',...
    'Units','centimeters',...
    'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*8.9 2*sqrt(2)*8.9]);
for ii = 1:length(Readout)            
    subplot(5,1,ii)
    hold on
    plot(x,y(ii,:))
    errorbar(x,yy(ii,:),yyerr(ii,:))
end
subplot(5,1,1)
title(['sin position'])
ylabel('std(predict(sin)-real(sin))')
xlabel('number of components')
subplot(5,1,2)
title('cos position')
ylabel('std(predict(cos)-real(cos))')
xlabel('number of components')
subplot(5,1,3)
title('Speed')
ylabel('std(predict(speed)-real(speed))')
xlabel('number of components')
subplot(5,1,4)
title('Predictor')
ylabel('% false predictions')
xlabel('number of components')
subplot(5,1,5)
title('position')
ylabel('std(predict(ang)-real(ang))')
xlabel('number of components')
%%
print(gcf, '-dpdf',num2str(kk))
kk = kk+1;

%% GPFA sectioned to rounds
% % % C = caim.S;C(caim.S_bin==0)= 0;
% % clear dat
% % roundstart = find(diff(scn.rounds)==1);
% % roundstart = [1;roundstart;size(C,2)];
% % if ~isempty(runstart)
% %     for i = 1:length(roundstart)-1
% %         dat(i).trialId = i;
% %         roundrun = false(1,size(C,2));
% %         roundrun(roundstart(i):roundstart(i+1)) = scn.running(roundstart(i):roundstart(i+1));
% %         dat(i).spikes = C(:,roundrun);
% %         dat(i).space = scn.distance(roundrun);
% %         dat(i).space(dat(i).space <1) = 1;
% %     end
% %     method = 'gpfa';
% %     xDim = 10;
% %     binWidth = 1; 
% % 
% %     % Extract neural trajectories
% %     result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
% %                         'kernSDList', kernSD,'binWidth',binWidth);
% %     [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
% % 
% %     figure('color',[1 1 1],...
% %         'renderer','painters',...
% %         'visible','off',...
% %         'Units','centimeters',...
% %         'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
% %         'PaperUnits','centimeters',...
% %         'PaperSize', [2*sqrt(2)*8.9 2*8.9])
% %     hold on
% %     mycolormap = hsv(max(space));
% % 
% %     for j = 1:size(seqTrain,2)
% %         aa = seqTrain(j).xorth';
% %         for i = 1:size(aa,1)-1
% %             plot3([aa(i,1) aa(i+1,1)],[aa(i,3) aa(i+1,3)],[aa(i,2) aa(i+1,2)],'color',mycolormap(round(dat(j).space(i)),:));         
% %         end
% %     end
% %     title('sectioned to rounds')
% %     hold off
% %     %%
% %     print(gcf, '-dpdf',num2str(kk))
% %     kk = kk+1;
   
%% sectioned to running periods
C = caim.S_bin;C(caim.S_bin==0)= 0;
C = smoother(C,kernSD,1);

clear dat
runstart = find(diff(scn.running)==1);
runstop = find(diff(scn.running)==-1);
if runstop(1)<runstart(1);runstop = runstop(2:end);end
if runstart(end)>runstop(end);runstart = runstart(1:end-1);end

method = 'pca';
xDim = round(size(C,1)/10);
if xDim >50; xDim = 50;end
if xDim <20; xDim = 20;end
% xDim = 30;
binWidth = 1; 
%     method = 'pca';
%     xDim = 10;%round((size(C,1)-1)/4);
%     kernSD = 5;
ii = true;
shrt = size(C,1);
while ii
    j = 1;
    for i = 1:length(runstart)
        int = runstart(i):runstop(i);
        if length(int)>20
            dat(j).trialId = j;
            dat(j).spikes = C(1:shrt,int);
            dat(j).space = scn.distance(int);
            dat(j).space(dat(j).space<1) = 1;
            j = j+1;
        end
    end

    % Extract neural trajectories
    result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                        'kernSDList', kernSD,'binWidth',binWidth);
    if ~isempty(result) || shrt <100
        ii = false;
    else
        shrt = shrt - 10;
    end
end
[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);

%%
if ~isempty(result)
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','off',...
        'Units','centimeters',...
        'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
        'PaperUnits','centimeters',...
        'PaperSize', [2*sqrt(2)*8.9 2*8.9])
    hold on
    mycolormap = hsv(max(space));

    for j = 1:size(seqTrain,2)
        aa = seqTrain(j).xorth';
    %     interp1(1:size(a,2),a(i,:),1:size(a,2)/size(caim.C,2):size(a,2),'spline');
        for i = 1:size(aa,1)-1
            plot3([aa(i,1) aa(i+1,1)],[aa(i,3) aa(i+1,3)],[aa(i,2) aa(i+1,2)],'color',mycolormap(round(dat(j).space(i)),:));         
        end
    end
    title('sectioned to running periods')
    hold off
    
    %%
    print(gcf, '-dpdf',num2str(kk))
    kk = kk+1;
    %% Cos distance to networks    
    PCAout = caim.PCAout;
    if isfield(caim.network,'netraster') && ~isempty(PCAout)
        dist = PCAout.sectioned.dist;
        distHist = PCAout.sectioned.distHist;
        confint = PCAout.sectioned.confint;
        distp = PCAout.sectioned.distp;
        
        clmplg = floor(256*[min(dist(:))/(min(dist(:))-max(dist(:))) 1-min(dist(:))/(min(dist(:))-max(dist(:)))]);
        mycolormap = zeros(256,3);
        j = clmplg(2);
        for i = 1:sum(clmplg)
            if i <clmplg(1)+2
                mycolormap(i,:) = [1 1-(clmplg(1)-(i-1))/clmplg(1) 1-(clmplg(1)-(i-1))/clmplg(1)];  
            else
                j = j-1;
                mycolormap(i,:) = [1-(clmplg(2)-j)/clmplg(2) 1-(clmplg(2)-j)/clmplg(2) 1];
            end
        end
        mycolormap(end,:) = [0 0 1];
        
        
        figure('color',[1 1 1],...
            'renderer','painters',...
            'visible','off',...
            'Units','centimeters',...
            'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
            'PaperUnits','centimeters',...
            'PaperSize', [2*8.9 2*sqrt(2)*8.9])
        
        x = .0:.005:1;
        subplot(2,1,1)
        plot(distHist(1,:),distHist(2,:))
        hold on
        plot(distHist(1,:),distHist(3,:))
        title(['Ks test: p = ' num2str(distp)])

        subplot(2,1,2)
        plotdist = dist;
        plotdist(dist > confint(1) & dist < confint(2)) = 0;
        plotdist(isnan(plotdist)) = 0;
        imagesc((plotdist)')
        colormap(mycolormap)
        colorbar
        %%
        print(gcf, '-dpdf',num2str(kk))
        kk = kk+1;
    end

end
    
%% centered around running onset
% %     clear dat
% %     win = -60:60;
% %     mycolormap = jet(length(win));
% %     runstart = find(diff(scn.running)==1);
% %     runstart = runstart(runstart>abs(win(1)) & runstart<size(C,2)-win(end));
% %     for i = 1:length(runstart)
% %         dat(i).trialId = i;
% %         dat(i).spikes = C(:,runstart(i)+win);
% %         dat(i).space = win;
% %     end
% %     
% %     method = 'gpfa';
% %     xDim = 10;
% %     binWidth = 1; 
% % 
% % 
% %     % Extract neural trajectories
% %     result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
% %                         'kernSDList', kernSD,'binWidth',binWidth);
% %     [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
% % 
% %     figure('color',[1 1 1],...
% %         'renderer','painters',...
% %         'visible','off',...
% %         'Units','centimeters',...
% %         'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
% %         'PaperUnits','centimeters',...
% %         'PaperSize', [2*sqrt(2)*8.9 2*8.9])
% %     hold on
% % 
% %     for j = 1:size(seqTrain,2)
% %         aa = seqTrain(j).xorth';
% %         for i = 1:size(aa,1)-1
% %             plot3([aa(i,1) aa(i+1,1)],[aa(i,3) aa(i+1,3)],[aa(i,2) aa(i+1,2)],'color',mycolormap(round(dat(j).space(i)+abs(win(1))+1),:));         
% %         end
% %     end
% %     title('sectioned around running onsets')
% %     hold off
% % end
% % %%
% % print(gcf, '-dpdf',num2str(kk))
% % kk = kk+1;

%% centered around network
C = caim.S;C(caim.S_bin==0)= 0;
clear dat
win = -35:35;
netpos = find(caim.network.netev);
netpos1 = netpos(netpos>abs(win(1)) & netpos<size(C,2)-win(end));

method = 'pca'; 
xDim = round(size(C,1)/10);
if xDim >30; xDim = 30;end
binWidth = 1; 

ii = true;
shrt = size(C,1);
while ii
    for i = 1:length(netpos1)
        dat(i).trialId = i;    
        dat(i).spikes = C(1:shrt,netpos1(i)+win);
        dat(i).space = win;
    end
    % Extract neural trajectories
    result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                        'kernSDList', kernSD,'binWidth',binWidth);
    if ~isempty(result) || shrt <100
        ii = false;
    else
        shrt = shrt - 10;
    end
end

[estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
%%
if ~isempty(result)
    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','off',...
        'Units','centimeters',...
        'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
        'PaperUnits','centimeters',...
        'PaperSize', [2*sqrt(2)*8.9 2*8.9])

    hold on
    win = 3:35;
    mycolormap = jet(length(win));
    maxplot =[1:22 24:40];%length(seqTrain);%30;
%     win = 1:69;
%     mycolormap = hsv(length(win));
%     mycolormap = [mycolormap(end/2:end,:);mycolormap(1:end/2,:)];
%     maxplot =[24:40 ];%length(seqTrain);%30;
    if size(seqTrain,2)<maxplot(end);maxplot = maxplot(1):size(seqTrain,2);end
    for j = maxplot
        aa = seqTrain(j).xorth';
        for i = win
            plot3([aa(i,1) aa(i+1,1)],[aa(i,2) aa(i+1,2)],[aa(i,3) aa(i+1,3)],'color',mycolormap(i-win(1)+1,:));         
        end
    end
    
    title('sectioned around networks')
    hold off
    %%
    print(gcf, '-dpdf',num2str(kk))
    kk = kk+1;
end

%% Cos distance network and running components   
PCAout = caim.PCAout;
if isfield(PCAout,'sectioned') && ~isempty(PCAout.sectioned.dual.distHist)
    dist = PCAout.sectioned.dual.dist;
    distHist = PCAout.sectioned.dual.distHist;
    confint = PCAout.sectioned.dual.confint;
    distp = PCAout.sectioned.dual.distp;

    clmplg = floor(256*[min(dist(:))/(min(dist(:))-max(dist(:))) 1-min(dist(:))/(min(dist(:))-max(dist(:)))]);
    mycolormap = zeros(256,3);
    j = clmplg(2);
    for i = 1:sum(clmplg)
        if i <clmplg(1)+2
            mycolormap(i,:) = [1 1-(clmplg(1)-(i-1))/clmplg(1) 1-(clmplg(1)-(i-1))/clmplg(1)];  
        else
            j = j-1;
            mycolormap(i,:) = [1-(clmplg(2)-j)/clmplg(2) 1-(clmplg(2)-j)/clmplg(2) 1];
        end
    end
    mycolormap(end,:) = [0 0 1];


    figure('color',[1 1 1],...
        'renderer','painters',...
        'visible','off',...
        'Units','centimeters',...
        'position',[10 1 [ 2*8.9 2*sqrt(2)*8.9]],...
        'PaperUnits','centimeters',...
        'PaperSize', [2*8.9 2*sqrt(2)*8.9])

    subplot(2,1,1)
    plot(distHist(1,:),distHist(2,:))
    hold on
    plot(distHist(1,:),distHist(3,:))
    title(['Ks test: p = ' num2str(distp)])

    subplot(2,1,2)
    plotdist = dist;
    plotdist(dist > confint(1) & dist < confint(2)) = 0;
    plotdist(isnan(plotdist)) = 0;
    imagesc((plotdist)')
    colormap(mycolormap)
    colorbar
    %%
    print(gcf, '-dpdf',num2str(kk))
    kk = kk+1;
end


%% centered around Airpuff
% if isfield(scn,'airpuff')
%     % C = caim.S;C(caim.S_bin==0)= 0;
%     stim = scn.airpuff;
%     clear dat
%     win = -30:120;
%     mycolormap = jet(length(win));
%     stimon = find(diff(stim.stimon)==1);
%     stimon = stimon(stimon>abs(win(1)) & stimon<size(C,2)-win(end));
%     
%     method = 'gpfa';
%     xDim = round(size(C,1)/10);
%     if xDim >30; xDim = 30;end
%     binWidth = 1; 
% 
%     ii = true;
%     shrt = size(C,1);
%     while ii
%         for i = 1:length(stimon)
%             dat(i).trialId = i;
%             dat(i).spikes = C(1:shrt,stimon(i)+win);
%             dat(i).space = win;
%         end
%         % Extract neural trajectories
%         result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
%                             'kernSDList', kernSD,'binWidth',binWidth);
%         if ~isempty(result) || shrt <100
%             ii = false;
%         else
%             shrt = shrt - 10;
%         end
%     end
%     [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
%     %%
%     if ~isempty(result)
%         figure('color',[1 1 1],...
%             'renderer','painters',...
%             'visible','off',...
%             'Units','centimeters',...
%             'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
%             'PaperUnits','centimeters',...
%             'PaperSize', [2*sqrt(2)*8.9 2*8.9])
%         hold on
%         win = 15:90;
%         mycolormap = jet(length(win));
%         maxplot = 30;
%         if length(size(seqTrain,2))<maxplot;maxplot = size(seqTrain,2);end
%         for j = 1:maxplot
%             aa = seqTrain(j).xorth';
%             for i = win%size(aa,1)-1
%                 plot3([aa(i,3) aa(i+1,3)],[aa(i,2) aa(i+1,2)],[aa(i,1) aa(i+1,1)],'color',mycolormap(i-win(1)+1,:));         
%             end
%         end
%         title('sectioned around airpuffs')
%         hold off
%         %%
%         print(gcf, '-dpdf',num2str(kk))
%         kk = kk+1;
%     end
% end
%%
if isfile([pathname filename '.pdf'])
    delete([pathname filename '.pdf'])
end
close all
files = dir('*.pdf');
append_pdfs([pathname filename '.pdf'],files.name)
% append_pdfs(['ICA CA1 example.pdf'],files.name)
delete(files.name)

end