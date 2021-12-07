%%
cclust = Airpuff.Fire.cclust;

meanf = 3;% 1 mean frequency of cell firing
% 2 mean frequency of cell firing during run
% 3 mean frequency of cell firing during rest
plcfld =  4;    % 4 placefield center
plcfldp = 5;    % 5 placefield p value
% 6 place coding vector angle
plcvct = 7;     % 7 place coding vector length
plcvctp = 8;    % 8 place coding vector p value
% 9 # response events
nstim = 10;% 10 # responded stimuli
% 11 percentage of network related firing
% 12 percentage of NON network related firing
netprob = 13;% 13 netprob
mpp = 14;% 14 mean MPP input

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

% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\correlations')
%%
% a = [];
a = cclust(:,meanf);
% a = cclust(cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05),meanf);
% a = cclust(cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05),meanf);
% a = cclust(cclust(:,nstim)<=2,meanf);
% a = cclust((cclust(:,nstim)<=2)|(cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05))|(cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05)),meanf);
% figure
% histogram(a)
% axis([0 30 0 50])
% hold on

b = cclust(cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05,meanf);
% figure
% histogram(b,30)
% axis([0 30 0 50])

c = cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,meanf);
% figure
% histogram(c,200)
% axis([0 30 0 50])
% hold off

d = cclust(cclust(:,nstim)>2,meanf);
% histogram(d,200)

boxplot([a' b' c' d'],[zeros(1,length(a)) ones(1,length(b)) 2*ones(1,length(c)) 3*ones(1,length(d))])

[p,t,stats] = kruskalwallis([a' b' c' d'],[zeros(1,length(a)) ones(1,length(b)) 2*ones(1,length(c)) 3*ones(1,length(d))]);
figure
cc = multcompare(stats);%,'CType','dunn-sidak');
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\histfire')
%%
% a = [];
a = cclust(:,netprob);
% a = cclust(cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05),netprob);
% a = cclust(cclust(:,plcvctp)>.05,netprob);
% histogram(a)
% hold on

b = cclust(cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05,netprob);
% histogram(b,30)
c = cclust(cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05,netprob);
% histogram(c,200)
hold off

d = cclust(cclust(:,nstim)>1,netprob);
% histogram(d,200)

boxplot([a' b' c' d'],[zeros(1,length(a)) ones(1,length(b)) 2*ones(1,length(c)) 3*ones(1,length(d))])

[p,t,stats] = kruskalwallis([a' b' c' d'],[zeros(1,length(a)) ones(1,length(b)) 2*ones(1,length(c)) 3*ones(1,length(d))]);

c = multcompare(stats);

%% Pie diagramm with place fields
figure('color',[ 1 1 1])
a = [];
if ~isnan(cclust(1,nstim)) %If stimulus was given
    % network, place, stim
    a(2) = length(find(cclust(:,netprob)>0 & (cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05) & cclust(:,nstim)>2))/length(cclust);
    % network, no place, stim
    a(1) = length(find(cclust(:,netprob)>0 & (cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05)) & cclust(:,nstim)>2))/length(cclust);
    % network, place, no stim
    a(3) = length(find(cclust(:,netprob)>0 & (cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05) & cclust(:,nstim)<=2))/length(cclust);
    % no network, place, stim
    a(4) = length(find(cclust(:,netprob)==0 & (cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05) & cclust(:,nstim)>2))/length(cclust);
    % no network, no place, stim
    a(5) = length(find(cclust(:,netprob)==0 & (cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05)) & cclust(:,nstim)>2))/length(cclust);
    % no network, place, no stim
    a(6) = length(find(cclust(:,netprob)==0 & (cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05) & cclust(:,nstim)<=2))/length(cclust);
    % network, no place, no stim
    a(7) = length(find(cclust(:,netprob)>0 & (cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05)) & cclust(:,nstim)<=2))/length(cclust);  
    % no network, no place, no stim
    a(8) = length(find(cclust(:,netprob)==0 & (cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05)) & cclust(:,nstim)<=2))/length(cclust);
    labels = {'net & plc & stim','net & stim','net & plc','plc & stim','stim','plc','net','non'};
else %if no stimulus was given
    % network, place
    a(1) = length(find(cclust(:,netprob)>0 & (cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05)))/length(cclust);
    % network, no place
    a(2) = length(find(cclust(:,netprob)>0 & (cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05))))/length(cclust);
    % no network, place
    a(3) = length(find(cclust(:,netprob)==0 & (cclust(:,plcfld)>0 & cclust(:,plcfldp)<=.05)))/length(cclust);
    % no network, no place
    a(4) = length(find(cclust(:,netprob)==0 & (cclust(:,plcfld)==0 | (cclust(:,plcfld)>0 & cclust(:,plcfldp)>.05))))/length(cclust);  
    labels = {'net & plc','net','plc','non'};
end
pie(a,labels)

% pie(a([1 2 3 5 6 7])/sum(a([1 2 3 5 6 7])),labels([1 2 3 5 6 7]))
% % exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\pieDombeck')

% bar([a; a],'stacked')
% ylim([0 .2])
% xlim([1.5 2.5])
% % exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\barLos')
%% Pie diagramm with place vector
figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]]) 
a = [];
if ~isnan(cclust(1,nstim)) %If stimulus was given
    % network, place, stim
    a(2) = length(find(cclust(:,netprob)>0 &  cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05 & cclust(:,nstim)>2))/length(cclust);
    % network, no place, stim
    a(1) = length(find(cclust(:,netprob)>0 &  (cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05)) & cclust(:,nstim)>2))/length(cclust);
    % network, place, no stim
    a(3) = length(find(cclust(:,netprob)>0 &  cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05 & cclust(:,nstim)<=2))/length(cclust);
    % no network, place, stim
    a(4) = length(find(cclust(:,netprob)==0 &  cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05 & cclust(:,nstim)>2))/length(cclust);
    % no network, no place, stim
    a(5) = length(find(cclust(:,netprob)==0 &  (cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05)) & cclust(:,nstim)>2))/length(cclust);
    % no network, place, no stim
    a(6) = length(find(cclust(:,netprob)==0 &  cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05 & cclust(:,nstim)<=2))/length(cclust);
    % network, no place, no stim
    a(7) = length(find(cclust(:,netprob)>0  &  (cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05)) & cclust(:,nstim)<=2))/length(cclust);
    % no network, no place, no stim
    a(8) = length(find(cclust(:,netprob)==0 &  (cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05)) & cclust(:,nstim)<=2))/length(cclust);
    labels = {'net & plc & stim','net & stim','net & plc','plc & stim','stim','plc','net','non'};
else %if no stimulus was given
    % network, place
    a(1) = length(find(cclust(:,netprob)>0 &  cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05 ))/length(cclust);
    % network, no place
    a(2) = length(find(cclust(:,netprob)>0 &   (cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05))))/length(cclust);
    % no network, place
    a(3) = length(find(cclust(:,netprob)==0 &  cclust(:,plcvct)>0 & cclust(:,plcvctp)<=.05 ))/length(cclust);
    % no network, no place
    a(4) = length(find(cclust(:,netprob)==0 &  (cclust(:,plcvct)==0 | (cclust(:,plcvct)>0 & cclust(:,plcvctp)>.05))))/length(cclust);  
    labels = {'net & plc','net','plc','non'};
end
b = a;
pie(a,labels)
% pie(a([1 2 3 5 6 7])/sum(a([1 2 3 5 6 7])),labels([1 2 3 5 6 7]))
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\pieLoso')