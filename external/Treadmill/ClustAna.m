function network = ClustAna(caim,scn)

network = caim.network;
if ~isfield(caim,'Y')
    return
end
if isempty(caim.network.netID) || size(caim.network.netraster,2)<4
    network.cellID = 0;
    network.netcorr = [];
    network.netcluster = [];
    network.clustConcept = [];
    network.clustbulk = [];
    network.clustfreq = [];  % Cluster frequency
    network.clustsize = [];
    network.cluststim = [];  % feature cells per cluster
    network.clustplace = [];
    network.clusttimes = [];
    network.clustcells = [];
    network.clustbulk = [];
    network.clustthresh = [];
    return
end
netraster = network.netraster;
netnum = network.netnum;

%% cluster sorting
% pearson's correlation
[netcorr,~] = corr(netraster');
a = ~isnan(diag(netcorr));
netcorr = netcorr(a,a);
Y = netcorr;
numcell = sum(a);

% % cos-dist
% celldist = zeros(size(netraster,1));
% for i = 1:length(celldist)
%     for j = 1:length(celldist)
%         celldist(i,j) = (netraster(j,:)*netraster(i,:)')./(sqrt(sum(netraster(j,:)))* sqrt(sum(netraster(i,:))) );
%     end
% end
% % numcell = size(netraster,1);
% % a = true(1,numcell);
% netcorr = celldist;
% Y = celldist(a,a);

Z = linkage(Y,'weighted','seuclidean');

%         D = pdist(Y,'seuclidean');     
%         c = optimalleaforder(Z,D);
[aa,bb] = histcounts((netnum),'normalization','probability');

sumaa = 0;
for ii = 1:length(aa)
    sumaa = sumaa + aa(ii);
    if sumaa >0.5
        meannet = bb(ii)+.5;
        break
    end
end 

numclust = round(numcell/meannet);

%% r-threshold for cluster number

% thresh = .3;
numit = 1000;
[~,thresh] = ShuffleClust(caim,numit,[],numclust,1);
%%
clustThresh = [];
k = 1;
num = -round(numclust/2):1:round(numcell)/2;

for j = num
    C = cluster(Z,'maxclust',numclust+j);
    B = zeros(1,max(C));

    for i = 1:max(C)        
        A = netcorr(C==i,C==i);
        A = tril(A,-1);
        A(A==0) = NaN;
        B(i) =  nanmean(A(:));
    end
    
    if nanmean(B) >thresh &&  k == 1
        newnum = numclust+j;
        clustsig = B;
        clustsig(2,:) = B>thresh;
        k = 0;
    end
   
    clustThresh(j+round(numclust/2)+1,:) = [numclust+j,nanmean(B)];   
end

if k ==1
    newnum = numclust+j;
    clustsig = B;
    clustsig(2,:) = B>thresh;
end
    
numclust = newnum;
%%
C = cluster(Z,'maxclust',numclust);
[~,c] = sort(C,'descend');
netcluster = Y(c,c);
cellID(a,1) = C;
cellID(a,2) = c;
cellID(~a,:) = NaN;

%% Shuffle analysis for 0-hypothethis of network size

[numshuf,~] = ShuffleClust(caim,numit,thresh,numclust,0);
clustp = sum(numshuf<=numclust)/numit;

%% concept cells in clusters

if isfield(scn,'airpuff')
    stim = scn.airpuff;
    isstim = false(1,length(a));
    isstim(stim.cellID(:,1)) = 1;
    isstim = isstim(a);
else
    isstim = false(1,length(a));
    isstim = isstim(a);
end

if isfield(scn,'cellID')
    isplace = false(1,length(a));
    isplace(scn.cellID(scn.cellID(:,6)<.05,1)) = 1;
    isplace = isplace(a);
else
    isplace = false(1,length(a));
    isplace = isplace(a);
end



if isfield(scn,'speedcorr') && ~isempty(scn.speedcorr)
    cellIDspeed = scn.speedcorr.cellID;
    isspeedpos = false(1,length(a));
    isspeedpos(cellIDspeed(:,1)>.9 & cellIDspeed(:,2)<.05 & cellIDspeed(:,3)<.05) = 1;
    isspeedpos = isspeedpos(a);
    isspeedneg = false(1,length(a));
    isspeedneg(cellIDspeed(:,1)>.9 & cellIDspeed(:,2)<.05 & cellIDspeed(:,3)<.05) = 1;
    isspeedneg = isspeedneg(a);
else
    isspeedpos = false(1,length(a));
    isspeedpos = isspeedpos(a);
    isspeedneg = false(1,length(a));
    isspeedneg = isspeedneg(a);
end

clustfreq = zeros(1,max(C));  % Cluster frequency
clustsize = zeros(1,max(C));
cluststim = zeros(1,max(C));  % feature cells per cluster
clustplace = zeros(1,max(C));
clustpos = zeros(1,max(C));
clustneg = zeros(1,max(C));
clusttimes = cell(1,max(C));
clustcells = cell(1,max(C));
clustbulk = cell(1,max(C));

bbb = (find(caim.network.netev ==1) >= (find(scn.network.stimon == 1,1,'first')))&(find(caim.network.netev ==1) <= (find(scn.network.stimon == 1,1,'last')));
%%
for i = 1:max(C)
    aa = C == i;
    clustsize(i) = sum(aa);
    clustfreq(i) = sum((sum(netraster(aa,:),1)>0))/size(netraster,2);
    cluststim(i) = sum(isstim(aa));
    clustplace(i) = sum(isplace(aa));
    clustpos(i) = sum(isspeedpos(aa));
    clustneg(i) = sum(isspeedneg(aa));
    clusttimes{i} = caim.network.netpos(sum(netraster(aa,:),1)>0);
    clustcells{i} = find(aa);
    if isfield(scn,'network') && isfield(scn.network,'bulkresp') && size(scn.network.bulkresp,1) == length(bbb)
        aaa = sum(netraster(aa,:),1)'>0 ;
        aaa = aaa(bbb);
        clustbulk{i} = scn.network.bulkresp(aaa,:);
    end 
end

%% Cluster analysis output
clustID.clustsig = clustsig;
clustID.isstim = cluststim ;
clustID.isplace = clustplace;
clustID.ispos = clustpos;
clustID.isneg = clustneg;

network.cellID = cellID;
network.netcorr = corr(netraster');
network.netcluster = netcluster;
network.clustThresh = clustThresh;
network.clustbulk = clustbulk;
network.clustfreq = clustfreq;  % Cluster frequency
network.clustsize = clustsize;
network.cluststim = cluststim;  % feature cells per cluster
network.clustplace = clustplace;
network.clustpos = clustpos;
network.clustneg = clustneg;
network.clusttimes = clusttimes;
network.clustcells = clustcells;
network.clustbulk = clustbulk;
network.clustthresh = thresh;
network.clustID = clustID;
network.numshuf = numshuf;
network.clustp = clustp;


%% cluster sorting with high active cells
[netcorr,~] = corr(netraster');
a = ~isnan(diag(netcorr)) & sum(netraster,2)>5;
netcorr = netcorr(a,a);
if isempty(netcorr) || length(netcorr)<3
    network.netHAcluster = [];return
end

Y = netcorr;
Z = linkage(Y,'weighted','seuclidean');
%         D = pdist(Y,'seuclidean');     
%         c = optimalleaforder(Z,D);
[aa,bb] = histcounts((netnum),'normalization','probability');
numcell = sum(a);
sumaa = 0;
for ii = 1:length(aa)
    sumaa = sumaa + aa(ii);
    if sumaa >0.5
        meannet = bb(ii)+.5;
        break
    end
end 
numclust = round(numcell/meannet);
if numclust<2;network.netHAcluster = []; return;end
C = cluster(Z,'maxclust',numclust);
[~,c] = sort(C,'descend');
netHAcluster = Y(c,c);

network.netHAcluster = netHAcluster;
return
%% network size change

cellID = true(size(netraster,1),1);
if isfield(scn,'cellID')
    cellID(scn.cellID(:,1)) = false;
end
if isfield(scn,'airpuff')
    stim = scn.airpuff;
    cellID(stim.cellID(:,1)) = false;
end
cellID = find(cellID);

netsize = NaN(length(cellID),size(netraster,2));
group = zeros(size(netsize));

for i = 1:length(cellID)
    a = netraster(cellID(i),:)>0;
    isstim = netraster(:,a);
    if sum(a)>3
        netsize(i,a) = sum(isstim,1)./max(sum(isstim,1));
        c = .05*round((1 : sum(a))*20/sum(a));
        group(i,a) = c;
    end
end
netnorm = zeros(sum(sum(~isnan(netsize))),2);
netnorm(:,1) = netsize(~isnan(netsize));
netnorm(:,2) = group(~isnan(netsize));
%%
if isfield(scn,'airpuff')
    stim = scn.airpuff;
    cellID = stim.cellID(:,1);
    netsize = NaN(length(cellID),size(netraster,2));
    group = zeros(size(netsize));
    for i = 1:length(cellID)        
        a = netraster(cellID(i),:)>0;
        isstim = netraster(:,a);
        if sum(a)>3
            netsize(i,a) = sum(isstim,1)./max(sum(isstim,1));
            c = .05*round((1 : sum(a))*20/sum(a));
            group(i,a) = c;
        end
    end
    netnormair = zeros(sum(sum(~isnan(netsize))),2);
    netnormair(:,1) = netsize(~isnan(netsize));
    netnormair(:,2) = group(~isnan(netsize));
else
    netnormair = [];
end
%%
if isfield(scn,'cellID')
    cellID = scn.cellID(:,1);
    netsize = NaN(length(cellID),size(netraster,2));
    group = zeros(size(netsize));
    for i = 1:length(cellID)
        a = netraster(cellID(i),:)>0;
        isstim = netraster(:,a);
        if sum(a)>3
            netsize(i,a) = sum(isstim,1)./max(sum(isstim,1));
            c = .05*round((1 : sum(a))*20/sum(a));
            group(i,a) = c;
        end
    end
    netnormplace = zeros(sum(sum(~isnan(netsize))),2);
    netnormplace(:,1) = netsize(~isnan(netsize));
    netnormplace(:,2) = group(~isnan(netsize));
else
    netnormplace = [];
end


%% network size output

network.netnorm = netnorm;
network.netnormair = netnormair;
network.netnormplace = netnormplace;



end