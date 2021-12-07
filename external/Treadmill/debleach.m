function aa = debleach(d,a)
%%
a = a-min(a);
a = a/max(a);

%% remoce artifacts if laser was instable
thresh = -5;
a = reshape(a,size(d));
c = zscore(a);  
if find(c<thresh)
    cc = d(c>thresh);
    c = c(c>thresh);
    c = interp1(cc,c,d);
    c(isnan(c))=mean(c(~isnan(c)));   
    a = c;
    
end

%% remove photobleaching

fs = length(d)/d(end);  
fc = [1]; % cutoff frequency 
Wn= fc/(fs/2);
n = 2;
[b,bb] = butter(n,Wn,'low');
aa = filtfilt(b,bb,a);
%%
ExpEqn = 'a*exp(-((x-b)/c))+d';
startPoints = [.2 0 1000 min(aa)];
upperBounds = [1 1 10000 .4];
lowerBounds =[.05 -1 500 0];
F1 = fit(d,aa,ExpEqn,'Start', startPoints,'upper',upperBounds,'lower',lowerBounds,'Exclude',d<100 & d>1000);
aa = F1(d);
aa = a-aa;
%%
% figure
% plot(d,a)
% hold on 
% plot(d,F1(d))
% plot(d,aa)
% hold off


end
