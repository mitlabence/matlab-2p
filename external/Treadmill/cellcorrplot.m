function cellcorrplot(cclust)
% cclust = Airpuff.Fire.cclust;

load('cclustID.mat','cclustID');
meanf = cclustID.meanf;             % mean frequency of cell firing
plcfld =  cclustID.plcfld;          % placefield center
% plcfldp = cclustID.plcfldp;       % placefield p value
plclength = cclustID.plclength;     % size of place field 
plcvang = cclustID.plcvctang;       % place coding vector angle
plcvct = cclustID.plcvct;           % place coding vector length
plcvctp = cclustID.plcvctp;         % place coding vector p value
nstim = cclustID.nstim;             % # responded stimuli
netperc = cclustID.netperc;         % percentage of network related firing
netprob = cclustID.netprob;         % netprob
mpp = cclustID.mpp;                 % mean MPP input
runperc = cclustID.runperc;         % percentage of running associated cells
runnerp = cclustID.runnerp;         % p value running
airpp = cclustID.airpp;             % p value number of response events
followcell = cclustID.followcell;   % logical if cell was in FOV of all experiments

%%
figure('position',[93 49 1722 948],'color',[1 1 1 ])
subplot(4,3,1)
scatter(cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,meanf),cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,plcvct))
[h,p] = corr(cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,meanf),cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,plcvct));
title(['Place-vector against activity, r = ' num2str(h) ', p = ' num2str(p)]);

subplot(4,3,2)
scatter(cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,mpp),cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,plcvct))
[h,p] = corr(cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05 & ~isnan(cclust(:,mpp)),mpp),cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05& ~isnan(cclust(:,mpp)),plcvct));
title(['Place-vector against Mpp, r = ' num2str(h) ', p = ' num2str(p)]);

subplot(4,3,3)
scatter(cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,nstim),cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,plcvct))
[h,p] = corr(cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,nstim),cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,plcvct));
title(['Place-vector against stim response, r = ' num2str(h) ', p = ' num2str(p)])

subplot(4,3,4)
scatter(cclust(:,meanf),cclust(:,netprob))
[h,p] = corr(cclust(:,meanf),cclust(:,netprob));
title(['Network against activity, r = ' num2str(h) ', p = ' num2str(p)]);

subplot(4,3,5)
scatter(cclust(:,mpp),cclust(:,netprob))
[h,p] = corr(cclust(~isnan(cclust(:,mpp)),mpp),cclust(~isnan(cclust(:,mpp)),netprob));
title(['Network against MPP, r = ' num2str(h) ', p = ' num2str(p)]);

if ~isnan(cclust(1,nstim))
    subplot(4,3,6)
    scatter(cclust(:,nstim),cclust(:,netprob))
    [h,p] = corr(cclust(:,nstim),cclust(:,netprob));
    title(['Network against stim response, r = ' num2str(h) ', p = ' num2str(p)])
    hold on
    aa = cclust(:,nstim);
    bb = cclust(:,netprob);
    bin = 100;
    binb = zeros(bin-1,1);
    bina = zeros(bin-1,1);
    bin = (min(aa):(max(aa)-min(aa))/bin:max(aa))';
    for i = 1:length(bin)-1
        bina(i) = mean([bin(i) bin(i+1)]);
        binb(i) = nanmean(bb(aa>=bin(i) & aa<=bin(i+1)));   
    end
    scatter(bina,binb,'linewidth',2)
    LnEqn = 'a*x+b';
    startPoints = [1 0];
    upperBounds = [10 2];
    % lowerBounds =[.05 -1 500 0];
    F1 = fit(bina(~isnan(binb)&bina<40&bina>0),binb(~isnan(binb)&bina<40&bina>0),LnEqn,'Start', startPoints,'upper',upperBounds);%,'lower',lowerBounds);
    h = get(gca,'ylim');
    plot(aa,F1(aa))
    set(gca,'ylim',h)
    hold off
end

subplot(4,3,7)
scatter(cclust(:,meanf),cclust(:,nstim))
[h,p] = corr(cclust(:,meanf),cclust(:,nstim));
title(['Stim response against activity, r = ' num2str(h) ', p = ' num2str(p)]);

subplot(4,3,8)
scatter(cclust(:,mpp),cclust(:,nstim))
[h,p] = corr(cclust(~isnan(cclust(:,mpp)),mpp),cclust(~isnan(cclust(:,mpp)),nstim));
title(['Stim response against MPP, r = ' num2str(h) ', p = ' num2str(p)]);

subplot(4,3,10)
scatter(cclust(:,meanf),cclust(:,mpp))
[h,p] = corr(cclust(~isnan(cclust(:,mpp)),meanf),cclust(~isnan(cclust(:,mpp)),mpp));
title(['MPP against activity, r = ' num2str(h) ', p = ' num2str(p)]);


subplot(4,3,11)
hold on
aa = cclust(cclust(:,mpp)>0,mpp);
bb = (1./(cclust(cclust(:,mpp)>0,meanf)));
aa(aa==inf)=nan;
bb(bb==inf)=nan;
bin = 100;
binb = zeros(bin-1,1);
bina = zeros(bin-1,1);
bin = (min(aa):(max(aa)-min(aa))/bin:max(aa))';
for i = 1:length(bin)-1
    bina(i) = mean([bin(i) bin(i+1)]);
    binb(i) = nanmean(bb(aa>=bin(i) & aa<=bin(i+1)));   
end
[h,p] = corr(bina(~isnan(binb)),binb(~isnan(binb)));
% [h,p] = corr(~isnan(aa),~isnan(bb));
scatter(aa,bb)
scatter(bina,binb,'linewidth',2)
title(['reciprocal activity against MPP, r = ' num2str(h) ', p = ' num2str(p)]);
LnEqn = 'a*x+b';
startPoints = [1 0];
upperBounds = [10 2];
% lowerBounds =[.05 -1 500 0];
F1 = fit(bina(~isnan(binb)&bina<3&bina>0),binb(~isnan(binb)&bina<3&bina>0),LnEqn,'Start', startPoints,'upper',upperBounds);%,'lower',lowerBounds);
h = get(gca,'ylim');
plot(aa,F1(aa))
set(gca,'ylim',h)
hold off

% exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\correlations')



