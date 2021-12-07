function caim = NetCorr(caim,scn)
if ~isfield(caim,'bulk') || ~isfield(caim,'Y')
    return
end

%%
fs = length(scn.tsscn)/scn.tsscn(end)*1000; 
Fs = 14.9;

d = scn.tsscn/1000;

a = sum(caim.C);
a = debleach(d,a);
a = a';
a = zscore(a);

b = sum(caim.S_bin);
b = b/10;

c = caim.bulk.trace(:,3);
c = zscore(c);

run = scn.running(1:end);

%%
netev = reshape(caim.network.netev,size(run));
bulkev = reshape(caim.bulk.bulkev,size(run));
caim.bulk.netfld(1,1) = length(find(netev(bulkev==1 & run ==0)))/length(find(netev(run ==0)));
caim.bulk.netfld(1,2) = length(find(netev(bulkev==0 & run ==0)))/length(find(netev(run ==0)));

%% calculate x-cross to get delay

win = 10;
alpha = .05;
max_lag = 5;

[ee,lags] = xcorr(a(run==0) , c(run==0),'coeff',round(win*fs));

% first line times, second line values
netxbulk = [(lags(lags/fs>-win & lags/fs<win))/fs;
            ee(lags/fs>-win & lags/fs<win)];
        
k = -round(win*fs)+1:round(win*fs)-1;
confint = 2./sqrt(sum(run==0)-abs(k));

[F,c_v,p] = granger_cause(a(run==0),c(run==0),alpha,max_lag);

netCp = NaN(2,2);
netCp(:,1) = [ee(lags==0)>confint(lags==0);c_v<F];

ee = max(ee(lags/fs>-(win-6) & lags/fs<(win-6)));
lags = lags(netxbulk(2,:)==ee); % max for shift in next step

% figure
% plot(netxbulk(1,:),netxbulk(2,:))
% hold on
% plot(netxbulk(1,:),confint)
% grid on

netxbulk = downSampResp(netxbulk);
%%

if sum(run==1) > 50 
    [ee1,lags1] = xcorr(a(run==1) , c(run==1),'coeff',round(win*fs));

    % first line times, second line values
    netxbulkrun = [(lags1(lags1/fs>-win & lags1/fs<win))/fs;
                ee1(lags1/fs>-win & lags1/fs<win)];
            
    confint = 2./sqrt(sum(run==1)-abs(k));
    [F1,c_v,p1] = granger_cause(a(run==1),c(run==1),alpha,max_lag);
    netCp(:,2) = [ee1(lags1==0)>confint(lags1==0);c_v<F1];
    ee1 = max(ee1(lags1/fs>-(win-6) & lags1/fs<(win-6)));
    lags1 = lags1(netxbulkrun(2,:)==ee1);
    netxbulkrun = downSampResp(netxbulkrun);
else
    netxbulkrun = nan(size(netxbulk));
    ee1 = NaN;
    lags1 = NaN;
    p1 = NaN;
    F1 = NaN;
end


%% calcualte pearson corr with corrected delay
% cc = zeros(1,length(c));
% if lags>=0
%     cc(1+lags:end) = c(1:end-lags);
% else
%     cc(1:end+lags) = c(1-lags:end);
% end
% [ee,p] = corr([a(run==0) ; cc(run==0)]');
%%

netcorr(1) = F;
netcorr(2) = p;
netcorr(3) = ee;
netcorr(4) = lags/fs;
netcorr(5) = F1;
netcorr(6) = p1;
netcorr(7) = ee1;
netcorr(8) = lags1/fs;

caim.bulk.netxbulk = netxbulk;
caim.bulk.netxbulkrun = netxbulkrun;
caim.bulk.netcorr = netcorr;
caim.bulk.netCp = netCp;

%%

alpha = .05;
max_lag = 5;
runs = find(diff(run));
if run(1) == 1, runs(1)= [];end
j = 1;
traces = nan(1,2*round(win*Fs)+1);
F = [];
c_v = [];
isrun = [];
issig = [];
for i = 1:length(runs)-1
    if length(runs(i):runs(i+1))>80     
        aa = a(runs(i):runs(i+1));
        cc = c(runs(i):runs(i+1));
        [ee,lags] = xcorr(aa,cc,'coeff',round(win*fs));
        k = round(lags*fs);
        confint = 2./sqrt(length(runs(i):runs(i+1)));
        [F(j),c_v(j)] = granger_cause(aa,cc,alpha,max_lag);
        issig(j,:) = [c_v(j)<F(j) max(ee(lags==0))>confint]; 
        isrun(j) = run(runs(i)+1)==1;       
        [temp] = downSampResp([lags;ee]);
        lags = temp(1,:);
        traces(j,:) = temp(2,:);
        j = j+1;
    end
end

% plot(lags/fs,mean(traces(~isrun & issig(:,2)' & issig(:,1)',:)))
% hold on
% plot(lags/fs,mean(traces(isrun & issig(:,2)' & issig(:,1)',:)))
% hold off

caim.bulk.nxbpart.traces = traces;
caim.bulk.nxbpart.lags = lags;
caim.bulk.nxbpart.isrun = isrun;
caim.bulk.nxbpart.issig = issig;
%% calculate individual x-cross

indcorr = zeros(size(caim.C,1),4);
for i = 1:size(caim.C,1)
    %%
    a = caim.C(i,:);
    a = zscore(a);
    [ee,lags] = xcorr(a(run==0) , c(run==0)); 
    lags = lags(max(ee)==ee);
    
    cc = zeros(1,length(c));
    if lags>=0
        cc(1+lags:end) = c(1:end-lags);
    else
        cc(1:end+lags) = c(1-lags:end);
    end
    [ee,p] = corr([a(run==0) ; cc(run==0)]');
    if p(2,1)<0.05 && ee(2,1)>0 && abs(lags)<300
        indcorr(i,1)=ee(2,1);
        indcorr(i,2)=p(2,1);
        indcorr(i,3)=lags/fs;
    end
end
indcorr(:,4) = caim.fireprop.fire(:,1);
caim.bulk.indcorr = indcorr;


function netxbulk = downSampResp(netxbulk)
%     Fs = 14.9;  
    time = netxbulk(1,:);
    if fs < 14
        timenew = time(1):.5*(1/fs):time(end);
        x = interp1(time,netxbulk(1,:),timenew,'spline');
        x(2,:) = interp1(time,netxbulk(2,:),timenew,'spline');   
        disp('X-correlogram is sampled up')
    elseif fs > 16
        timenew = time(1):2*(1/fs):time(end); 
        x = interp1(time,netxbulk(1,:),timenew,'spline');
        x(2,:) = interp1(time,netxbulk(2,:),timenew,'spline');
        disp('X-correlogram is sampled down')
    else
        x = netxbulk;
    end
    setlgth = 2*win*Fs+1;
    if length(x) > setlgth
        fact = (length(x)-1)/2-(setlgth-1)/2;
    else
        fact = 0;
    end
    netxbulk = x(:,1+fact:end-fact);
end

end