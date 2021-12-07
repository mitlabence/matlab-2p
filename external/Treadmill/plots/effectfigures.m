
mouseID = cell2mat(Airpuff.Response.airpuff.mouseID(:,2));
% trial = 3;
effects = Airpuff.Response.airpuff.effects;
runtime = Airpuff.Response.airpuff.runtime;
waittime = Airpuff.Response.airpuff.waittime;
distance = [];
waittimetemp = [];
resp = [];
nonresp = [];
pupil = [];
bulk = [];
aporder = [];
for i = 1:size(effects,1)
    if mouseID(i) == trial
        a = cell2mat(effects(i));
        aporder = [aporder; (1:size(a,1))'];
        pupiltemp = a(:,4,2)-mean(a(:,4,1));
        pupiltemp = pupiltemp/max(pupiltemp);
        pupil =  [pupil; pupiltemp];
        bulktemp = a(:,3,2)-mean(a(:,3,1));
        bulktemp = bulktemp/max(bulktemp);
        bulktemp(bulktemp==0) = NaN;
        bulk = [bulk; bulktemp];
        nonresptemp = (a(:,2,2)-mean(a(:,2,1)));
        nonresp = [nonresp; nonresptemp];
        resptemp = a(:,1,2)-mean(a(:,1,1));%+a(:,1,3);%-a(:,1,1);
        resp = [resp; resptemp];   
        b = cell2mat(runtime(i));
        b(:,2) =b(:,2)./max(b(:,2));
        distance = [distance; b];
        waittimetemp = [waittimetemp; cell2mat(waittime(i))];
    end
end

a =  ~isnan(distance(:,1)) & ~isnan(waittimetemp)  & aporder<100 ;
distance = distance(a,2);
aporder = aporder(a);
resp = resp(a);
nonresp = nonresp(a);
pupil = pupil(a);
bulk = bulk(a);
waittimetemp = waittimetemp(a);
% scatter(aporder,nonresp,'g')
% hold on
% scatter(aporder,runtimetemp,'b')
% scatter(aporder,bulk,'r')
% hold off


%%

x = [aporder resp nonresp distance waittimetemp pupil bulk];%
[r,p] = corr(x);
r(r==1) = NaN;
r(p>.05) = NaN;
% r = abs(r);
% b = imagesc(r);
% set(b,'AlphaData',~isnan(r))
% colorbar
% ax = gca;
% ax.XTickLabel = {'APorder' 'resp' 'nonresp' 'distance' 'waittime' 'pupil' 'bulk'};
% ax.YTickLabel = ax.XTickLabel;
%%
aa = resp;
bin = 100;round(max(aporder)/3);
bina = zeros(bin-1,1);
binaporder = zeros(bin-1,1);
binresp = zeros(bin-1,1);
binnonresp = zeros(bin-1,1);
bindistance = zeros(bin-1,1);
binpupil = zeros(bin-1,1);
binbulk = zeros(bin-1,1);
bin = (min(aa):(max(aa)-min(aa))/bin:max(aa))';
for i = 1:length(bin)-1
    bina(i) = nanmean([bin(i) bin(i+1)]);
    binaporder(i) = nanmean(aporder(aa>=bin(i) & aa<=bin(i+1)));
    binresp(i) = nanmean(resp(aa>=bin(i) & aa<=bin(i+1)));
    binnonresp(i) = nanmean(nonresp(aa>=bin(i) & aa<=bin(i+1)));
    bindistance(i) = nanmean(distance(aa>=bin(i) & aa<=bin(i+1)));
    binpupil(i) = nanmean(pupil(aa>=bin(i) & aa<=bin(i+1)));
    binbulk(i) = nanmean(bulk(aa>=bin(i) & aa<=bin(i+1)));
end



[rr(1),pp(1)] = corr(bina(~isnan(binresp)),binresp(~isnan(binresp)));
[rr(2),pp(2)] = corr(bina(~isnan(binnonresp)),binnonresp(~isnan(binnonresp)));
[rr(3),pp(3)] = corr(bina(~isnan(binbulk)),binbulk(~isnan(binbulk)));
[rr(4),pp(4)] = corr(bina(~isnan(bindistance)),bindistance(~isnan(bindistance)));
[rr(5),pp(5)] = corr(bina(~isnan(binpupil)),binpupil(~isnan(binpupil)));


% subplot(5,
% scatter(bina,binresp,'g','linewidth',1)
% hold on
% scatter(bina,binnonresp,'MarkerEdgeColor',[0 .5 0])
% scatter(bina,binbulk,'r')
% scatter(bina,bindistance,'b')
% scatter(bina,binpupil,'y')
% hold off
% [rr,pp] = corr([binaporder binresp binnonresp  bindistance binpupil binbulk]);
%%
mouseID = cell2mat(Airpuff.Response.airpuff.mouseID(:,2));
effects = Airpuff.Response.airpuff.effects;
runtime = Airpuff.Response.airpuff.runtime;
waittime = Airpuff.Response.airpuff.waittime;
distance = cell(3,1);
waittimetemp = cell(3,1);
resp = [];
pupil = [];
for trial = 1:3
    for i = 4:size(effects,1)
        if mouseID(i) == trial
            a = cell2mat(effects(i));
            a = a(:,1,2);%+a(:,1,3);%-a(:,1,1);
            resp = [resp; a];   
            b = cell2mat(runtime(i));
            b = nanmean(b(1:20,2));
            % b(:,2) =b(:,2)./max(b(:,2));
            distance{trial} = [distance{trial}; b];
            b = cell2mat(waittime(i));
            b= nanmean(b);
            waittimetemp{trial} = [waittimetemp{trial}; b];
        end
    end
end

a = [distance{1} distance{2} distance{3}];
% a = [waittimetemp{1} waittimetemp{2} waittimetemp{3}];
plot(a','color',[.5 .5 .5])
hold on
plot(mean(a),'b')
hold off
%% Pooled Data
lnwd = 1;
ftsz = 3;
stim = Airpuff.Response.airpuff;
win = 1:121;
sumresp = zeros(size(stim.sumresp));
nresp = zeros(size(stim.nresp));
trial = 1;
for i = 1:size(stim.nums,1)
    sumresp(i,:) = stim.plcresp(i,:)/stim.nums(i,1);%/stim.nums(i,2);
    nresp(i,:) = stim.spcresp(i,:)/stim.nums(i,1);%/stim.nums(i,2);
end

% sumresp = sumresp(cell2mat(stim.mouseID(:,2))==3,:);

exclude = [ ];
for i = 1:length(exclude)
    temp = find(cell2mat(stim.mouseID(:,4))==exclude(i));
    for j = 1:length(temp)
        stim.mouseID{temp(j),2} = 0;
    end
end
excludebulk = [ ];
for i = 1:length(excludebulk)
    temp = find(cell2mat(stim.mouseIDbulk(:,4))==excludebulk(i));
    for j = 1:length(temp)
        stim.mouseIDbulk{temp(j),2} = 0;
    end
end

x = stim.times(1,win)/1000;

subplot(7,1,1)


y = nanmean(sumresp(cell2mat(stim.mouseID(:,2))==trial,win),1);
yy = nanmean(nresp(cell2mat(stim.mouseID(:,2))==trial,win),1);
plot(x,y,'color',[0 1 0],'linewidth',lnwd)
hold on
b = std(sumresp(cell2mat(stim.mouseID(:,2))==trial,win),1)./sqrt(size(sumresp(cell2mat(stim.mouseID(:,2))==trial,win),1));
fill([x,fliplr(x)],[y-b,fliplr(y+b)],[0 .7 0],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)

plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd)
axis tight
axis off
hold off

subplot(7,1,2)
plot(x,yy,'color',[0 .7 0],'linewidth',lnwd)
hold on

b = nanstd(nresp(cell2mat(stim.mouseID(:,2))==trial,win),1)./sqrt(size(nresp(cell2mat(stim.mouseID(:,2))==trial,win),1));
fill([x,fliplr(x)],[yy-b,fliplr(yy+b)],[0 .7 0],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)

plot([0 0],[min(yy-b) max(yy+b)],'m','linewidth',lnwd)
axis tight
axis off
hold off

subplot(7,1,3)
y = nanmean(stim.bulk(cell2mat(stim.mouseIDbulk(:,2))==trial,win),1);
b = nanstd(stim.bulk(cell2mat(stim.mouseIDbulk(:,2))==trial,win),1)./sqrt(size(stim.bulk(cell2mat(stim.mouseIDbulk(:,2))==trial,win),1));
plot(x,y,'r','linewidth',lnwd)
hold on
fill([x,fliplr(x)],[y-b,fliplr(y+b)],[.7 0 0],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)
plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd)
axis tight
axis off
hold off

subplot(7,1,4)
y = stim.speed(cell2mat(stim.mouseID(:,2))==trial,win);
b = std(y,1)./sqrt(size(y,1));
y = mean(y,1);
plot(x,y,'color',[0 176/255,240/255],'linewidth',lnwd)
hold on
fill([x,fliplr(x)],[y-b,fliplr(y+b)],[0 176/255,240/255],...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.3)
plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd)
axis tight
axis off
hold off

subplot(7,1,5)
y = stim.pupil(cell2mat(stim.mouseID(:,2))==trial,win);
b = std(y,1)./sqrt(size(y,1));
y = mean(y,1);
plot(x,y,'y','linewidth',lnwd)
hold on
fill([x,fliplr(x)],[y-b,fliplr(y+b)],'y',...
    'EdgeColor',[1 1 1],...
    'EdgeAlpha',0,...
    'FaceAlpha',.2)
plot([0 0],[min(y-b) max(y+b)],'m','linewidth',lnwd)
axis tight
axis off
hold off
% 
% title('a)','Units','centimeter','position',[-1 6.5],'FontSize',ftsz)
%%
% y = c;
% % y = (~isnan(y));
% x = 1:size(y,2);
% yerr = nanstd(y)/sqrt(size(y,1));
% y = nanmean(y);
% 
% LnEqn = 'a*x+b';
% startPoints = [1 y(1)];
% % upperBounds = [1 1 10000 .4];
% % lowerBounds =[.05 -1 500 0];
% F1 = fit(x',y',LnEqn,'Start', startPoints,'Weights',1./yerr');%,'upper',upperBounds,'lower',lowerBounds);
% 
% figure
% hold on
% errorbar(x,y,yerr,'*')
% plot(x,F1(x))
% hold off

%% summary figure

% figure('color',[1 1 1],...
%     'position',[300 100 2*[594 420]],...
%     'visible','on',...
%     'PaperUnits','centimeters',...
%     'Units','centimeters',...
%     'renderer','painters')
% 
% for ii = 1:3
% trial = ii;
% effectfigures
% subplot(5,3,ii)
% scatter(bina,binresp,'g','linewidth',1)
% title(['Trial ' num2str(trial)])
% ax = gca;
% ax.XLim = [0 max(bina)];
% text(ax.XLim(2)*2/3,ax.YLim(2),['r = ' num2str(round(rr(1),2)) ', p = ' num2str(round(pp(1),2))])
% if ii == 1;ylabel('\Delta resp');end
% 
% subplot(5,3,ii+3)
% scatter(bina,binnonresp,'MarkerEdgeColor',[0 .5 0])
% ax = gca;
% ax.XLim = [0 max(bina)];
% text(ax.XLim(2)*2/3,ax.YLim(2),['r = ' num2str(round(rr(2),2)) ', p = ' num2str(round(pp(2),2))])
% if ii == 1;ylabel('\Delta non-resp');end
% 
% subplot(5,3,ii+6)
% scatter(bina,binbulk,'r')
% ax = gca;
% ax.XLim = [0 max(bina)];
% text(ax.XLim(2)*2/3,ax.YLim(2),['r = ' num2str(round(rr(3),2)) ', p = ' num2str(round(pp(3),2))])
% if ii == 1;ylabel('Bulk Amp');end
% 
% subplot(5,3,ii+9)
% scatter(bina,bindistance,'b')
% ax = gca;
% ax.XLim = [0 max(bina)];
% text(ax.XLim(2)*2/3,ax.YLim(2),['r = ' num2str(round(rr(4),2)) ', p = ' num2str(round(pp(4),2))])
% if ii == 1;ylabel('Run dist');end
% 
% subplot(5,3,ii+12)
% scatter(bina,binpupil,'y')
% ax = gca;
% ax.XLim = [0 max(bina)];
% text(ax.XLim(2)*2/3,ax.YLim(2),['r = ' num2str(round(rr(5),2)) ', p = ' num2str(round(pp(5),2))])
% if ii == 1;ylabel('Pupil');end
% xlabel('AP');
% end
% 
% % printpdf('C:\Users\Admin\Dropbox\Dentate in-vivo Project\Corr to DG resp',1)