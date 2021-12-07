function aa = dfoverf(d,x)
%%
a = x;
a = a-min(a);
a = a/max(a);

%% remoce artifacts if laser was instable
thresh = -3;
a = reshape(a,size(d));
c = zscore(a);  
if find(c<thresh)
    cc = d(c>thresh);
    c = a(c>thresh);
    c = interp1(cc,c,d);
    c(isnan(c))=mean(c(~isnan(c)));   
    c = c-min(c);
    c = c/max(c);
    a = c;
end

%% define f0 as the low pass filtered signal

fs = length(d)/d(end);  
fc = [.01]; % cutoff frequency 1
Wn= fc/(fs/2);
n = 2;
[b,bb] = butter(n,Wn,'low');
aa = filtfilt(b,bb,a);

aa = a./aa;

%%

% figure
% plot(d,aa)

end