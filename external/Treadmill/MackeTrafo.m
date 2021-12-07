function [Simil,SimShuf] = MackeTrafo(caim,scn)

E = logical(caim.network.netev);
Z = logical(caim.S_bin);
running = scn.running;
X = caim.C./caim.Df;
% X = caim.S;X(caim.S_bin==0)= 0;
if sum(E)<5
    Simil = [];
    SimShuf = [];
    return
end
%% running related activity
runstart = find(diff(scn.running)==1);
runstop = find(diff(scn.running)==-1);

if runstop(1)<runstart(1);runstop = runstop(2:end);end
if runstart(end)>runstop(end);runstart = runstart(1:end-1);end

clear dat

j = 1;
for i = 1:length(runstart)
    int = runstart(i):runstop(i);
    if length(int)>0
        dat(j).trialId = j;
        dat(j).spikes = X(:,int);
        dat(j).space = scn.distance(int);
        dat(j).space(dat(j).space<1) = 1;
        j = j+1;
    end
end
datrun = dat;
Xrun = [dat.spikes];
%%
runIdx = 2;
binWidth = 1; 
kernSD = 5;
method = 'pca';
xDim = size(X,1);
result = neuralTraj(runIdx, datrun, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD,'binWidth',binWidth);
                
[estParams, ~] = postprocess(result, 'kernSD', kernSD);
NetScore = estParams.Corth;
hasSpikesBool1 = result.hasSpikesBool;
sRun = zeros(1,size(X,1));
sRun(hasSpikesBool1) = estParams.explained;


hasSpikes1 = find(hasSpikesBool1);
Newscore = zeros(xDim);
for i = 1:size(NetScore,2)
    Newscore(hasSpikesBool1,hasSpikes1(i)) =  NetScore(:,i);
end
hasSpikesBool1(:) = true;
Yrun = Newscore;

%% network related activity

win = -30:30;
netpos = find(E);
netpos = netpos(netpos>abs(win(1)) & netpos<size(X,2)-win(end));

clear dat
for i = 1:length(netpos)
    dat(i).trialId = i;    
    dat(i).spikes = X(:,netpos(i)+win);
    dat(i).space = win;
end
datnet = dat;
Xnet = [dat.spikes];
%%
runIdx = 2;
binWidth = 1; 
kernSD = 5;
method = 'pca';
xDim = size(X,1);
result = neuralTraj(runIdx, datnet, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD,'binWidth',binWidth);
                
[estParams, ~] = postprocess(result, 'kernSD', kernSD);
NetScore = estParams.Corth;
hasSpikesBool1 = result.hasSpikesBool;
sNet = zeros(1,size(X,1));
sNet(hasSpikesBool1) = estParams.explained;


hasSpikes1 = find(hasSpikesBool1);
Newscore = zeros(xDim);
for i = 1:size(NetScore,2)
    Newscore(hasSpikesBool1,hasSpikes1(i)) =  NetScore(:,i);
end
hasSpikesBool1(:) = true;
Ynet = Newscore;

%% Similarity factors (Krzanowsky 1979)
L = Ynet;
M = Yrun;
if size(L,2)>size(M,2)
    L = L(:,1:size(M,2));
elseif size(L,2)<size(M,2)
    M = M(:,1:size(L,2));
end
%     S = L*M'*M*L';
% S = (L*M')^2;
S = zeros(size(L,2));
for j = 1:size(M,2)
    a = L(:,j);
    for i = 1:size(L,2) 
        b = M(:,i);
        if sqrt(sum(a.^2))*sqrt(sum(b.^2)) ~= 0
            S(i,j) = ((a'*b)./(sqrt(sum(a.^2))*sqrt(sum(b.^2))))^2;
        end
    end
end
    
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

Simil =  [sum(diag(S)) sum(Er)^2];   

%%
numit = 1000;
SimShuf = zeros(numit,2);

parfor ii = 1:numit
    datrand = randomX(X,E,running);
    %%
    runIdx = 2;
    binWidth = 1; 
    kernSD = 5;
    method = 'pca';
    xDim = size(X,1);
    result = neuralTraj(runIdx, datrand, 'method', method, 'xDim', xDim,... 
                        'kernSDList', kernSD,'binWidth',binWidth);

    [estParams, ~] = postprocess(result, 'kernSD', kernSD);
    NetScore = estParams.Corth;
    hasSpikesBool1 = result.hasSpikesBool;
    sRand = zeros(1,size(X,1));
    sRand(hasSpikesBool1) = estParams.explained;


    hasSpikes1 = find(hasSpikesBool1);
    Newscore = zeros(xDim);
    for i = 1:size(NetScore,2)
        Newscore(hasSpikesBool1,hasSpikes1(i)) =  NetScore(:,i);
    end
    hasSpikesBool1(:) = true;
    Yrand = Newscore;

    %% Similarity factors (Krzanowsky 1979)
    L = Yrand(:,:);
    M = Yrun(:,:);
    if size(L,2)>size(M,2)
        L = L(:,1:size(M,2));
    end
%     S = L*M'*M*L';
    S = zeros(size(L,2));
    for j = 1:size(M,2)
        a = L(:,j);
        for i = 1:size(L,2) 
            b = M(:,i);
            if sqrt(sum(a.^2))*sqrt(sum(b.^2)) ~= 0
                S(i,j) = ((a'*b)./(sqrt(sum(a.^2))*sqrt(sum(b.^2))))^2;
            end
        end
    end
    
    %% Eros Approach (Yang 2004)
    w = zeros(size(sRand));    
    for i = 1:length(w)
        w(i) = mean([sRand(i) sRun(i)]);        
    end
    w = w./sum(w);
    Er = zeros(size(sRand));
    for i = 1:length(w)
        Er(i) = w(i)*abs(L(:,i)'*M(:,i));
    end
    
    SimShuf(ii,:) =  [sum(diag(S)) sum(Er)^2];   
end

end

function [datrand,Xrand] = randomX(X,E,running)
    rest = find(running==0);
    
    %% random network points during rest
    Erand = false(size(E));
    Erand(randsample(rest,sum(E))) = true;
%     Erand = E;
    %% random reshuffle of data
%     Xrest = X(:,rest);
%     for i = 1:size(X,1)
%         tpt = randperm(size(Xrest,2),1);
%         Xrest(i,:) = Xrest(i,[tpt:end 1:tpt-1]); 
%     end
%     X(:,rest) = Xrest;
    %% create shuffled dataset
    win = -30:30;
    netpos = find(Erand);
    netpos = netpos(netpos>abs(win(1)) & netpos<size(X,2)-win(end));

    for i = 1:length(netpos) 
        dat(i).trialId = i;
        dat(i).spikes = X(:,netpos(i)+win);
        dat(i).space = win;
    end
    datrand = dat;
    Xrand = [dat.spikes];
end

%%
% save('M226.mat','E','running','Z','X','Xrun','Xnet','Xrand')