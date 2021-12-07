function PCAout = PCAana(caim,scn,method)
if ~isfield(caim,'S')
    PCAout = [];
    return
end
addpath(genpath('/media/2Photon/Matlab code/Martin Source Code/gpfa_v0203'))
% C = caim.C; C = C./caim.Df;     
C = caim.S;C(caim.S_bin==0)= 0;
clear dat
dat.trialId = 1;
dat.spikes = C;
%%
% [a,score,~,~,explained] = pca(C);
%%
% Select method to extract neural trajectories:
% 'gpfa' -- Gaussian-process factor analysis
% 'fa'   -- Smooth and factor analysis
% 'ppca' -- Smooth and probabilistic principal components analysis
% 'pca'  -- Smooth and principal components analysis
if strcmp(method,'pca')
    disp('Good old PCA')
    runIdx = 2;
    method = 'pca';
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
    if sum(hasSpikesBool)<xDim
        xDim = sum(hasSpikesBool);
    end
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
    kernSD = 5;
    xDim = 30;
    yOut = smoother(C, kernSD,1);
    b = rica(yOut',xDim,'IterationLimit',1000);
    % b = rica(a,q,'IterationLimit',1000);
    a = yOut'*b.TransformWeights;
    score = b.TransformWeights;
    hasSpikesBool = true(size(C,1),1);
end

%% Components triggered to running
% scn.totdist(run),1:10,abs(a(run,1:10)

%% Components triggered to network
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

if isfield(caim.network,'netraster')
    netraster = caim.network.netraster;
    netraster = netraster(hasSpikesBool,:);
    [dist,shufdist,confint,distHist,distPerc] = PCAshuf(netraster,hasSpikesBool,score,xDim);
    
    a1 = abs(dist(~isnan(dist(:))));
    a2 = abs(shufdist(~isnan(shufdist(:))));
    [~,p]  = kstest2(a1,a2);
    
    PCAout.dist = dist;
    PCAout.distp = p;
    PCAout.confint = confint;
    PCAout.distHist = distHist;
    PCAout.distPerc = distPerc;
else
    PCAout.dist = [];
    PCAout.distp = NaN; 
    PCAout.confint = [];
    PCAout.distHist = [];
    PCAout.distPerc = [];
end
%%
explained = var(a);

linEqn = 'a+b*(x)';
int = 4:round(length(explained)/3);
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


% a1 = zeros(size(a));
% a2 = a1;
% for i = 1:size(abinrun,2)
%     a1(:,i) = smooth(a(:,i).*abs(abinrun(:,i)),1);
%     a2(:,i) = smooth(a(:,i).*abs(abinnet(:,i)),1);
% end


a1 = zeros(size(a));
a2 = a1;
for i = 1:size(abinrun,2)
%     a1(:,i) = smooth(a(:,i).*abinrun(:,i),500);
%     a2(:,i) = smooth(abs(a(:,i).*abinnet(:,i)),500);
    a1(:,i) = cumsum(abs(abinrun(:,i)))/sum(abs(abinrun(:,i)));
    a2(:,i) = cumsum(abs(abinnet(:,i)))/sum(abs(abinnet(:,i)));
end


PCAout.abinrun = sum(abs(abinrun),2);
PCAout.abinnet = sum(abs(abinnet),2);
PCAout.cumrun = nanmean(a1,2);
PCAout.cumnet = nanmean(a2,2);
[PCAout.delay,PCAout.delaybin] = histcounts(adist(:,:,:)/15,[-300:1:300],'normalization','probability');

%% Sectioned PCA
if 0%strcmp(method,'pca') && sum(scn.running)>150 
    PCAout.sectioned = PCAsectioned(caim,scn,'pca');
    [Simil,SimShuf] = MackeTrafo(caim,scn);
    PCAout.sectioned.Simil  = Simil;
    PCAout.sectioned.SimShuf = SimShuf;
else
    PCAout.sectioned = [];
    PCAout.sectioned.Simil  = [];
    PCAout.sectioned.SimShuf = [];
end

end