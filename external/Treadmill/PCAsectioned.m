function out = PCAsectioned(caim,scn,method)
if ~isfield(caim,'C')
    out = [];
    return
end
% C = caim.C; C = C./caim.Df;     
C = caim.S;C(caim.S_bin==0)= 0;
if strcmp(method,'gpfa')
%     C = smoother(C,5,1);
    xDim = round(size(C,1)/10);
    if xDim >50; xDim = 50;end
    if xDim <20; xDim = 20;end
elseif strcmp(method,'pca')
    xDim = round(size(C,1))-1;
end

%% sectioned to running periods
    

runstart = find(diff(scn.running)==1);
runstop = find(diff(scn.running)==-1);
if isempty(runstart)
    out = [];
    return
end

if runstop(1)<runstart(1);runstop = runstop(2:end);end
if runstart(end)>runstop(end);runstart = runstart(1:end-1);end


binWidth = 1; 
runIdx = 2;
kernSD = 5;
ii = true;
shrt = size(C,1);
clear dat
dat.trialId = [];
dat.spikes = [];
dat.space = [];

while ii
    j = 1;
    for i = 1:length(runstart)
        int = runstart(i):runstop(i);
        if length(int)>0
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
    if ~isempty(result)
        ii = false;
    elseif shrt <=100
        break
    else
        shrt = shrt -10;
    end
    
%     if ~isempty(result)
%         ii = false;
%     elseif xDim <=3
%         break
%     else
%         xDim = xDim-1;
%     end
    
end
%%

if isfield(caim.network,'netraster') && ~isempty(result)
    [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
    score = estParams.Corth;
    hasSpikesBool = result.hasSpikesBool;
    netraster = caim.network.netraster;
    sRun = zeros(1,size(C,1));
    sRun(hasSpikesBool) = estParams.explained;
    
    testallcells = 1;
    if testallcells == 1
        hasSpikes = find(hasSpikesBool);
        Newscore = zeros(size(netraster,1),size(netraster,1));
        for i = 1:size(score,2)
            Newscore(hasSpikesBool,hasSpikes(i)) =  score(:,i);
        end
        score = Newscore;
        hasSpikesBool(:) = true;
    else
        netraster = netraster(hasSpikesBool,:);
    end
    
    [dist,shufdist,confint,distHist,distPerc,pInd] = PCAshuf(netraster,hasSpikesBool,score,xDim);  
    a1 = abs(dist(~isnan(dist(:))));
    a2 = abs(shufdist(~isnan(shufdist(:))));
    [~,p]  = kstest2(a1,a2);
    [rankp]  = ranksum(a1,a2);
    %%    
    out.dist = dist;
    out.confint = confint;
    out.score = score;
    out.distp = p;
    out.rankp = rankp;
    out.hasSpikesBool = hasSpikesBool;
    out.distHist = distHist;
    out.distPerc = distPerc;
    out.pInd = pInd;
else
    out.dist = [];
    out.confint = NaN;
    out.score = [];
    out.distp = [];
    out.hasSpikesBool = [];
    out.distHist = [];
    out.distPerc = [];
    out.pInd = [];
end

%% centered around network
clear dat
win = -30:30;
netpos = find(caim.network.netev);
netpos1 = netpos(netpos>abs(win(1)) & netpos<size(C,2)-win(end));
if length(netpos1) >5
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

    %%
    [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);
    NetScore = estParams.Corth;
    hasSpikesBool1 = result.hasSpikesBool;
    sNet = zeros(1,size(C,1));
    sNet(hasSpikesBool1) = estParams.explained;


    hasSpikes1 = find(hasSpikesBool1);
    Newscore = zeros(size(score));
    for i = 1:size(NetScore,2)
        Newscore(hasSpikesBool1,hasSpikes1(i)) =  NetScore(:,i);
    end
    hasSpikesBool1(:) = true;
    NetScore = Newscore;
    %% shuffle analysis
    [dist,shufdist,confint,distHist,distPerc,pInd] = PCAshuf(NetScore,hasSpikesBool1,score,xDim);  
    a1 = abs(dist(~isnan(dist(:))));
    a2 = abs(shufdist(~isnan(shufdist(:))));
    [~,p2]  = kstest2(a1,a2);
    [rankp]  = ranksum(a1,a2);
%     distHist = histcounts(a1,0:.01:1,'normalization','probability');
%     distHist(2,:) = histcounts(a2,0:.01:1,'normalization','probability');
    %% Similarity factors (Krzanowsky 1979)
    L = NetScore;
    M = score;
    if size(L,2)>size(M,2)
        L = L(:,1:size(M,2));
    end
    S = L*M'*M*L';
    Simil =  sum(diag(S));
    %% Eros Approach (Yang 2004)
    w = zeros(size(sNet));    
    for i = 1:length(w)
        w(i) = mean([sNet(i) sRun(i)]);        
    end
    w = w./sum(w);
    Er = zeros(size(sNet));
    for i = 1:length(w)
        Er(i) = w(i)*abs(L(:,i)'*M(:,i));
    end
    Simil(1,2) = sum(Er);
    %%
    out.dual.dist = dist;
    out.dual.distp = p2;
    out.dual.rankp = rankp;
    out.dual.RunScore = score;
    out.dual.NetScore = NetScore;
    out.dual.hasSpikesBool = hasSpikesBool;
    out.dual.confint = confint;
    out.dual.distHist = distHist;
    out.dual.distPerc = distPerc;
    out.dual.pInd = pInd;
    out.dual.simil = Simil;
else
    out.dual.dist = [];
    out.dual.distp = NaN;
    out.dual.rankp = []; 
    out.dual.RunScore = score;
    out.dual.NetScore = [];
    out.dual.hasSpikesBool = [];
    out.dual.confint = [];
    out.dual.distHist = [];
    out.dual.distPerc = [];
    out.dual.pInd = [];
    out.dual.simil = [];
end
    