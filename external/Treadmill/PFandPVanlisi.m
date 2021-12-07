numice = [1:11];
% The rows refer to experimental conditions the were pooled. Decide here
% which experiments you want to pool
numexp = [5  6];

% Since we have no idea how big the files will be in each trial, we keep
% the variable empty here and youse grouping variables for later analysis.
pv = [];
pvWT = [];
pvKA = [];
pvWA = [];
pvICA = [];
pf = [];
pfWT = [];
pfKA = [];
pfWA = [];
pfICA = [];
groupExp = [];
load('cclustID.mat')
% This loop goes through the rows (the experiemts)
for i = 1:length(numexp)
    
    % first sort out mice that might not contain the desired variable or
    % are excluded
    fullCAIM = numice;
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if ~isfield(CAIM(numexp(i),fullCAIM(j)),'cclust')
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];    
    
    % Here happens the read out
    % The loop goes through the included mice
    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);          
        cclust = CAIM(numexp(i),k).cclust;
        pf = [pf; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plclength)];
        pv = [pv; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)];
        
        if i == 1 && k<6
            pfWT = [pfWT; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plclength)];
            pvWT = [pvWT; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)];
            groupExp = [groupExp; 1*ones(length(cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)),1)];
        elseif i == 1 && k>=6
            pfKA = [pfKA; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plclength)];
            pvKA = [pvKA; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)];
            groupExp = [groupExp; 2*ones(length(cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)),1)];
        elseif i == 2 && k>=6
            pfICA = [pfICA; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plclength)];
            pvICA = [pvICA; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)];
            groupExp = [groupExp; 3*ones(length(cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)),1)];
        elseif i == 2 && k<6
            pfWA = [pfWA; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plclength)];
            pvWA = [pvWA; cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)];
            groupExp = [groupExp; 4*ones(length(cclust(cclust(:,cclustID.plcvctp)<.05,cclustID.plcvct)),1)];
        end
    end
    
end

 figure('color',[1 1 1],...
            'renderer','painters',...
            'visible','on',...
            'Units','centimeters',...
            'position',[10 1 [ 2*8.9 2*8.9]],...
            'PaperUnits','centimeters',...
            'PaperSize', [2*8.9 2*8.9]);
        
subplot(2,2,1)
boxplot(pv,groupExp)
ax = gca;
ax.XTickLabel = {'WT' 'KA' 'KA - ICA' 'WT - ICA'};
ylabel('PV length')
title('Place tuning vector length')

subplot(2,2,2)
boxplot(pf,groupExp)
ax = gca;
ax.XTickLabel = {'WT' 'KA' 'KA - ICA' 'WT - ICA'};
ylabel('PF length')
title('Place field length')

subplot(2,2,3)
x = 0:.1:1.9;
y = cumsum(histcounts(pvWT(~isnan(pvWT)),x,'normalization','probability'));
plot(x(2:end),y)
hold on
y = cumsum(histcounts(pvKA(~isnan(pvKA)),x,'normalization','probability'));
plot(x(2:end),y)
y = cumsum(histcounts(pvICA(~isnan(pvICA)),x,'normalization','probability'));
plot(x(2:end),y)
y = cumsum(histcounts(pvWA(~isnan(pvWA)),x,'normalization','probability'));
plot(x(2:end),y)
hold off
xlabel('PV length')
ylabel('Cumulativ Probability')
legend({'WT' 'KA' 'KA - ICA' 'WT - ICA'})

subplot(2,2,4)
x = 10:10:400;
y = cumsum(histcounts(pfWT,x,'normalization','probability'));
plot(x(1:end-1),y)
hold on
y = cumsum(histcounts(pfKA,x,'normalization','probability'));
plot(x(1:end-1),y)
y = cumsum(histcounts(pfICA,x,'normalization','probability'));
plot(x(1:end-1),y)
y = cumsum(histcounts(pfWA,x,'normalization','probability'));
plot(x(1:end-1),y)
hold off
xlabel('PF length')
ylabel('Cumulativ Probability')
legend({'WT' 'KA' 'KA - ICA' 'WT - ICA'})
%%
[p,tbl,stats]  = anova1(pv,groupExp);
figure
multcompare(stats);
%%
[p,tbl,stats]  = anova1(pf,groupExp);
figure
multcompare(stats);