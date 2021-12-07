%% CellIDplot
Cues = load('Z:\Martin Pofahl\Cues.mat');
Baseline = load('Z:\Martin Pofahl\Baseline.mat');
Airpuff = load('Z:\Martin Pofahl\Airpuff.mat');
load('Z:\Martin Pofahl\BigFatCluster.mat');
load('cclustID.mat')


%% frequencies base
num = 3:5;
firefreq = [];
netfreq = [];
netoccur = [];
% shuffreq = [];
shufoccur = [];

for i = num%1:size(CAIM,1)   
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
        
        firefreq = [firefreq; CAIM(i,k).fireprop.meanfire(2,1:3)];
        netoccur = [netoccur; CAIM(i,k).network.netfreq(1,1:3)];   
        netfreq = [netfreq; CAIM(i,k).network.netfreq(2,1:3)];
        shufoccur = [shufoccur; CAIM(i,k).network.ShufMean(1,1:3)];
    end  
end

n = 9;
fire = [mean(firefreq);
    std(firefreq)/sqrt(n)];
network = [mean(netfreq);
    std(netfreq)/sqrt(n)];
%% frequencies cue
num = 8;
firefreq = [];
netfreq = [];
netoccur = [];
% shuffreq = [];
shufoccur = [];

for i = num%1:size(CAIM,1)   
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
        
        firefreq = [firefreq; CAIM(i,k).fireprop.meanfire(2,1:3)];
        netoccur = [netoccur; CAIM(i,k).network.netfreq(1,1:3)];   
        netfreq = [netfreq; CAIM(i,k).network.netfreq(2,1:3)];
        shufoccur = [shufoccur; CAIM(i,k).network.ShufMean(1,1:3)];
    end  
end

n = 9;
fire = [mean(firefreq);
    std(firefreq)/sqrt(n)];
network = [mean(netfreq);
    std(netfreq)/sqrt(n)];

%% cell IDs

use = [3 4 5 8];% 
mice = [103 155 158 194 195 224 226 227 234];% Mice to use cells from
% mice = [103 158 194 195 226 227 234];
cellID = cclust(:,cclustID.expID,:);
cellID = max(cellID,[],3);
cellID = cellID == mice;
cellID = max(cellID,[],2);


a = zeros(length(use),8);
b = zeros(length(use),1);
n = zeros(length(use),1);

for i = 1:length(use)
    if use(i)<=5
        n(i) = sum(cell2mat(Baseline.Fire.mouseID(:,2))==use(i));
    elseif use(i)<=8
        n(i) = sum(cell2mat(Cues.Fire.mouseID(:,2))==use(i)-5);
    end
    
    cclusttemp = cclust(cellID,:,use(i));
    cclusttemp = cclusttemp(cclusttemp(:,cclustID.expID)>0,:);
    isnet = cclusttemp(:,cclustID.netprob)>0;
    isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
    isspeed = cclusttemp(:,cclustID.speedcorr)>.9 & cclusttemp(:,cclustID.speedcorrp2)<=.05;
    % network, place, no speed
    a(i,1) = sum(isnet &  isplace &  ~isspeed);
    % no network, place, no speed
    a(i,2) = sum(~isnet &  isplace &  ~isspeed);
    % network, place, speed 
    a(i,3) = sum(isnet &  isplace &  isspeed);
    % no network, place, speed
    a(i,4) = sum(~isnet &  isplace &  isspeed);
    % network, no place, speed 
    a(i,5) = sum(isnet & ~isplace &  isspeed);
    % no network, no place, speed
    a(i,6) = sum(~isnet &  ~isplace &  isspeed);
    % network, no place, no speed
    a(i,7) = sum(isnet  &  ~isplace & ~isspeed);
    % no network, no place, no speed
    a(i,8) = sum(~isnet &  ~isplace & ~isspeed);
    b(i) = length(cclusttemp);   
end
n(:) = 1;


aa = [sum(a(1:3,:),1); a(4,:)];
bb = [sum(b(1:3)); b(4)];
aa = 100*(aa./bb);

% place cells
PC = sum(aa(:,1:2),2);
PC(:,2) = aa(:,1)./PC;
% place % speed cells
PSC = sum(aa(:,3:4),2);
PSC(:,2) = aa(:,3)./PSC;
% speed cells
SC = sum(aa(:,5:6),2);
SC(:,2) = aa(:,5)./SC;
