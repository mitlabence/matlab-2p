% How to extract pooled data from the big CAIM variable 

% The coloums of CAIM refer to individual mice that were pooled. Set the following value to
% decide which mice should be included:
numice = [1:4];
% The rows refer to experimental conditions the were pooled. Decide here
% which experiments you want to pool
numexp = [5];

%% Read out skalar variable

% Define the variable where youwant to store the data:
numrounds = nan(length(numexp),size(CAIM,2)); 
% Allocate the space with nan so that you can easily ingore non existing
% entries in later analysis.
% In this example we pool the number of ran rounds. You can pool any
% variable!

% This loop goes through the rows (the experiemts)
for i = 1:length(numexp)
 
    % first sort out mice that might not contain the desired variable or
    % are excluded
    fullCAIM = numice;
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        % This if-statement checks whether the variable you want to extract exists in th trials.
        % Otherwise it is skipt
        if isempty(CAIM(numexp(i),fullCAIM(j)).behave)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];    
    
    % Here happens the read out
    % The loop goes through the included mice
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        numrounds(i,k) = CAIM(numexp(i),k).behave.numrounds; 
    end
end

% No you can do plotting, for analysis, copying to excel (what ever) with
% your data
figure
subplot(1,2,1)
plot(numrounds','color',[.5 .5 .5])
hold on
boxplot(numrounds)
xlabel('mouse ID')
ylabel('# rounds')
subplot(1,2,2)
plot(numrounds,'color',[.5 .5 .5])
hold on
boxplot(numrounds')
xlabel('Session')
ylabel('# rounds')

%% Read out 2D variable

% Define the variable where youwant to store the data:
placefields = cell(length(numexp),size(CAIM,2));
% If your data has more dimensions and they are different for every mouse,
% you can use a cell variable for pooling.

% This loop goes through the rows (the experiemts)
for i = 1:length(numexp)
    % first sort out mice that might not contain the desired variable or
    % are excluded
    fullCAIM = numice;
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if ~isfield(CAIM(numexp(i),fullCAIM(j)),'plcfield') || isempty(CAIM(numexp(i),fullCAIM(j)).plcfield)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];    
    
    % Here happens the read out
    % The loop goes through the included mice
    figure('position',[100 300 1700 500])
    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        placefieldTemp = CAIM(numexp(i),k).plcfield;
        placefields(i,k) = {placefieldTemp};
        
        % Inside of this loop you can do also some analysis or plotting.
        subplot(1,size(CAIM,2),k)
        imagesc(placefieldTemp)       
    end
    
    colormap(jet)
end

%% Accumulate data from individual cells

% Define the variable where youwant to store the data:
meanf = [];
% Since we have no idea how big the files will be in each trial, we keep
% the variable empty here and youse grouping variables for later analysis.
groupExp = [];
groupMou = [];

% This loop goes through the rows (the experiemts)
for i = 1:length(numexp)
    
    % first sort out mice that might not contain the desired variable or
    % are excluded
    fullCAIM = numice;
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if ~isfield(CAIM(numexp(i),fullCAIM(j)),'fireprop') || isempty(CAIM(numexp(i),fullCAIM(j)).fireprop)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];    
    
    % Here happens the read out
    % The loop goes through the included mice
    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        
        meanfTemp = CAIM(numexp(i),k).fireprop.fire(:,1);
        meanf = [meanf; meanfTemp];
        groupExp = [groupExp; i*ones(length(meanfTemp),1)];
        groupMou = [groupMou; k*ones(length(meanfTemp),1)];
        
       
    end
    
end

figure
subplot(1,2,1)
boxplot(meanf,groupMou)
xlabel('mouse ID')
ylabel('# rounds')
subplot(1,2,2)
boxplot(meanf,groupExp)

%% read from cclust
load('cclustID.mat')
% Define the variable where youwant to store the data:
CellCode = [];
plcvctp = [];
plcinfop = [];
% This loop goes through the rows (the experiemts)
for i = 1:length(numexp)
    
    % first sort out mice that might not contain the desired variable or
    % are excluded
    fullCAIM = numice;
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if ~isfield(CAIM(numexp(i),fullCAIM(j)),'cclust') || isempty(CAIM(numexp(i),fullCAIM(j)).cclust)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];    
    
    % Here happens the read out
    % The loop goes through the included mice
    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        
        cclust = CAIM(numexp(i),k).cclust;
        isplace = cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05;
        isspeed = cclust(:,cclustID.speedcorr)>.9 & cclust(:,cclustID.speedcorrp2)<=.05;
        CellCode(i,1,k) = sum(isplace &  ~isspeed)./size(cclust,1);
        CellCode(i,2,k) = sum(isplace &  isspeed)./size(cclust,1);
        CellCode(i,3,k) = sum(~isplace & isspeed)./size(cclust,1);
        CellCode(i,4,k) = sum(~isplace & ~isspeed)./size(cclust,1);
        plcvctp = [plcvctp; cclust(:,cclustID.plcvctp)];
        plcinfop = [plcinfop; cclust(:,cclustID.spatialInfop)];       
    end
    
end
plcvctp(isnan(plcvctp)) = [];
plcinfop (isnan(plcinfop)) = [];
sigperc = [sum(plcvctp<=.05 & plcinfop<=.05) sum(plcvctp>.05 & plcinfop<=.05) sum(plcvctp<=.05 & plcinfop>.05) sum(plcvctp>.05 & plcinfop>.05)];
figure
subplot(1,2,1)
scatter((plcinfop),(plcvctp))
xlabel('p value skaggs')
ylabel('p value place-vector')
title('p-value correlation')

subplot(1,2,2)
scatter((plcinfop),(plcvctp))
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('p value skaggs')
ylabel('p value place-vector')
title('p-value correlation log scale')

% print(gcf, '-dpdf', 'Skaggs info correlation'); 
%%
figure
bar(permute(CellCode,[3 2 1]),'stacked')
legend({'place','place & speed','speed','none'},'location','eastoutside')
legend('boxoff')