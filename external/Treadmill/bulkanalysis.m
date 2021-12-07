function caim = bulkanalysis(caim,scn)
if isfield(caim,'bulk') && ~isempty(caim.bulk)
    bulk = caim.bulk;
else
    return
end
%%
fs = length(scn.tsscn)/scn.tsscn(end)*1000;              % Sampling frequency   
if isfield(caim.bulk,'traceMEC')
    disp('MEC input signal found')
    x = caim.bulk.traceMEC(1,:);
elseif isfield(caim.bulk,'traceMC')
    disp('Mossy cell input signal found')
    x = caim.bulk.traceMC(1,:);
else
    disp('No bulk signal found')
    return
end
 
% x = debleach(scn.tsscn/1000,x);
% 
% if length(find(scn.running==1))>100
%     x = (x -mean(x(scn.running==0)))/std(x(scn.running==0));
% else
%     x = zscore(x);
% end
x = dfoverf(scn.tsscn/1000,x);

%% low pass filter

n = 2;
fc = [.1]; % cutoff frequency
Wn= fc/(fs/2);
[b,a] = butter(n,Wn,'low');
lpx = filtfilt(b,a,x);

% figure
% win = 1:5000;
% plot(scn.tsscn(win)/1000,x(win))
% hold on
% plot(scn.tsscn(win)/1000,.2*scn.running(win))
% plot(scn.tsscn(win)/1000,lpx(win))
% hold off

%% high pass filter

n = 2;
fc = [.01]; % cutoff frequency
Wn= fc/(fs/2);
[b,a] = butter(n,Wn,'high');
hpx = filtfilt(b,a,x);

% figure
% hold on
% % plot(scn.tsscn(win)/1000,x(win))
% plot(scn.tsscn(win)/1000,2*scn.running(win))
% plot(scn.tsscn(win)/1000,hpx(win))
%% bandpass filter

n = 2;
fc = [.1 2]; % cutoff frequency
Wn= fc/(fs/2);
[b,a] = butter(n,Wn,'bandpass');
bpx = filtfilt(b,a,x);

% figure
% win = 1:5000;
% plot(scn.tsscn(win)/1000,scn.running(win))
% hold on
% plot(scn.tsscn(win)/1000,bpx(win))
% hold off

%%

if length(find(scn.running==1))>100
    [fftspec,ratio,baseline,speedcorr,bulkev,bincorr] = runana;
else
    [fftspec,ratio,baseline,speedcorr,bulkev,bincorr] = standana;
end

bulk.trace = [x hpx bpx lpx];
bulk.baseline = baseline;
bulk.fftspec = fftspec;
bulk.ratio = ratio;
bulk.speedcorr = speedcorr;
bulk.bincorr = bincorr;
bulk.bulkev = bulkev;

caim.bulk = bulk;
%%

function [fftspec,ratio,baseline,speedcorr,bulkev,bincorr] =runana
        
    baseline = [mean(x) mean(x(scn.running==1)) mean(x(scn.running==0));
        var(x)  var(x(scn.running==01)) var(x(scn.running==0))];
    
    %% Correlation speed versus amplitude
    
    speed = scn.speed*100;
    speed(speed<-20) = 0;
    speed(end+1)=speed(end);
    speed = smooth(speed,30);
    
    win = 5;
    
    [ee,lags] = xcorr(zscore(x(scn.running==1)) , (speed(scn.running==1)),'coeff',round(win*fs));
    
    lags = lags(max(ee(lags/fs>-100 & lags/fs<100))==ee);
    
    speedlag = zeros(1,length(speed));
    if lags>=0
        speedlag(1+lags:end) = speed(1:end-lags);
    else
        speedlag(1:end+lags) = speed(1-lags:end);
    end
    
    [ee,p] = corr([hpx(scn.running==1)'  ;  speedlag(scn.running==1)]');
%     [ee,p] = corr([hpx(scn.running==1)'  ;  speed(scn.running==1)']');
    
    % correlation coefficient
    speedcorr(1) = ee(2,1); 
    % p value
    speedcorr(2) = p(2,1);
    % time lag
    speedcorr(3) = lags/fs;
    
    speed = speed(1:length(x));
    speed(speed<=0) = 0;
    bin = min(speed) :max(speed)/20: max(speed);
    a = zeros(1,length(bin)-1);
    aerr = zeros(1,length(bin)-1);
    for i = 1:length(bin)-1
        a(i) = mean(x(speed>bin(i)&speed<=bin(i+1)));
        aerr(i) = std(x(speed>bin(i)&speed<bin(i+1)))/sqrt(sum(speed>bin(i)&speed<bin(i+1)));
    end
    
    bincorr = [bin(1:end-1); a; aerr];
    
%     scatter(speed(scn.running==1),x(scn.running==1))
%     hold on
%     errorbar(bin(1:end-1),a,aerr)
%     hold off
    %% Fourier Transform to identify peaks

    z = hpx;                      
    L = length(z);        % Length of signal

    Z = fft(z);
    P = abs(Z/L);
    P = P(1:floor(L/2+1));
    P(2:end-1) = 2*P(2:end-1);
    f = fs*(0:(L/2))/L;

    z1 = hpx(scn.running==0);               
    L = length(z1);        % Length of signal

    Z = fft(z1);
    P1 = abs(Z/L);
    P1 = P1(1:floor(L/2+1));
    P1(2:end-1) = 2*P1(2:end-1);
    f1 = fs*(0:(L/2))/L;

    z2 = hpx(scn.running==1);                      
    L = length(z2);        % Length of signal

    Z = fft(z2);
    P2 = abs(Z/L);
    P2 = P2(1:floor(L/2+1));
    P2(2:end-1) = 2*P2(2:end-1);
    f2 = fs*(0:(L/2))/L;  

    % figure
    % plot(f,smooth(P,100)) 
    % if ~isempty(f2)
    %     hold on
    %     plot(f1,smooth(P1,100)) 
    %     plot(f2,smooth(P2,100))  
    %     hold off
    % end
    
    %% Spectral ratio of runnning versus not running 
    
    [P,Q] = rat(length(P2)/length(P1));
    ynew = smooth(resample(P1,P,Q),400);
    ynew2 = smooth(P2,400);
    if length(ynew) > length(ynew2)
        ynew = ynew(1:length(ynew2));
    elseif length(ynew) < length(ynew2)
        ynew2 = ynew(1:length(ynew));
    end
    ratio(1) = min(abs(ynew2(f2>.1 & f2<4)./ynew(f2>.1 & f2<4)));
    if fs > 10
        ratio(2) = max(abs(ynew2(f2>4 & f2<15)./ynew(f2>4 & f2<15)));
    else 
        ratio(2) = 0;
    end
    fftspec = [reshape(f2,size(ynew))  ynew ynew2];
    
%     figure
%     hold on
%     plot(f2,ynew)   
%     plot(f2,ynew2)
%     plot(f2(f2>.1 & f2<4),abs(ynew2(f2>.1 & f2<4)./ynew(f2>.1 & f2<4)))
%     plot(f2(f2>4 & f2<15),abs(ynew2(f2>4 & f2<15)./ynew(f2>4 & f2<15)))
%     plot(f2,abs(ynew2./ynew),...
%         'linewidth',2,...
%         'color',[1 0 0])
%     set(gcf,'color',[0 0 0])
%     set(gca,'Xlim',[0 15],...
%             'fontsize',20,...
%             'YColor',[1 1 1],...
%             'XColor',[1 1 1],...
%             'Color',[0 0 0],...
%             'linewidth',2)
%      grid on
%      plot([0, 15],[1, 1],'--w')
%     
%      hold off
%      exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\bulk fft ratio')

    %% Try to make binary trace

    runon = diff(scn.running);
    runoff = find(runon==-1);
    runon = find(runon==1);

    if runon(1)<runoff(1) 
        runon = runon(2:end);
    end
    if length(runon) < length(runoff)
        runoff = runoff(1:end-1);
    elseif length(runon) < length(runoff)
        runoff = runoff(2:end);
    end
    
    bulkev = zeros(1,length(x));
    for i = 1:length(runoff)
        %%
        xx = bpx(runoff(i):runon(i));
        yy = zeros(length(xx),1);
        thresh = 0;
        yy(xx>thresh) = 1;
        bulkev(runoff(i):runon(i)) = yy;

    end
    
end

function [fftspec,ratio,baseline,speedcorr,bulkev,bincorr]  = standana
    baseline = [mean(x) mean(x) 0;
        std(x)  std(x) 0];
    fftspec = [];
    ratio = [0 0];
    speedcorr = [NaN NaN NaN];
    bincorr = NaN(3,20);
    xx = bpx;
    yy = zeros(length(xx),1);
    thresh = 0;
%     yy(xx<thresh) = thresh;
    yy(xx>thresh) = 1;
    bulkev = yy;
end

end