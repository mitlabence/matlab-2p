%%
stim = Airpuff.Response.airpuff;
Effects = stim.effects;
trial = 4;
Effects = Effects(cell2mat(stim.mouseID(:,2))==trial);
% Effects = Effects(4:end);
Effects = Effects(1:end-1);
lgth = 18;
a = zeros(length(Effects),lgth);
b = zeros(length(Effects),lgth);
c = zeros(length(Effects),lgth);
d = zeros(length(Effects),lgth);

figure
hold on
for i = 1:length(Effects)%[1:3 6:9]%%
    effects = Effects{i};
    effects(:,1,1) = mean(effects(:,1,1));
    effects(:,2,1) = mean(effects(:,2,1));
    effects = effects(:,:,2)-effects(:,:,1);%effects(:,:,3)+
    a(i,:) = effects(1:lgth,1);
    b(i,:) = effects(1:lgth,2);
    c(i,:) = effects(1:lgth,3);
    d(i,:) = effects(1:lgth,4);
%     plot(effects(:,1),'g')
%     plot(effects(:,2),'b')
%     plot(effects(:,3),'r')
%     plot(effects(:,4),'y')
    
end

c(c==0) = nan;
plot(nanmean(a),'g')

plot(nanmean(b),'b')
plot(nanmean(c),'r')
plot(mean(d),'y')
hold off

exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\Fascilitaion')
%%
y = c;
% y = (~isnan(y));
x = 1:size(y,2);
yerr = nanstd(y)/sqrt(size(y,1));
y = nanmean(y);

LnEqn = 'a*x+b';
startPoints = [1 y(1)];
% upperBounds = [1 1 10000 .4];
% lowerBounds =[.05 -1 500 0];
F1 = fit(x',y',LnEqn,'Start', startPoints,'Weights',1./yerr');%,'upper',upperBounds,'lower',lowerBounds);

figure
hold on
errorbar(x,y,yerr,'*')
plot(x,F1(x))
hold off

%% binned correlation plot  
figure
hold on
cc = 1:9;%
% cc = [6 8 9];%
aa = reshape(a(cc,:),[length(cc)*lgth 1]);
bb = reshape(b(cc,:),[length(cc)*lgth 1]);

bin = 50;
binb = zeros(bin-1,1);
bina = zeros(bin-1,1);
bin = (min(aa):(max(aa)-min(aa))/bin:max(aa))';
for i = 1:length(bin)-1
    bina(i) = mean([bin(i) bin(i+1)]);
    binb(i) = nanmean(bb(aa>=bin(i) & aa<=bin(i+1)));
    
end
[h,p] = corr(bina(~isnan(binb)),binb(~isnan(binb)));
% [h,p] = corr(aa,bb);
title(['h = ' num2str(h) ', p = ' num2str(p)]);
scatter(aa,bb)
scatter(bina,binb,'linewidth',2)
hold off