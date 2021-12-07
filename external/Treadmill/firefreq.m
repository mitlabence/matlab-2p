function fireprop = firefreq(caim,scn)
%%

if isstruct(caim) && isfield(caim,'S_bin')
    S = caim.S_bin;
    C = caim.C;
    Df = caim.Df;
else
    fireprop = [];
    return
end

tsscn = scn.tsscn;
run = scn.running;

Fs = 1000*length(tsscn)/tsscn(end);            % Sampling frequency
T = 1/Fs;

lges   = tsscn(end)/1000/60;                    %total time in min
lrun   = T*length(run(run==1))/60;              %time running
lstand = T*length(run(run==0))/60;              %time standing

fire = zeros(size(S,1),3);
firetot = zeros(size(S,1),3);
meanfire = zeros(3,3);
meanamplitude = zeros(2,3);
amplitude = zeros(size(S,1),3);


%% 
    
fire(:,1) = sum(S,2)/lges;
fire(:,2) = sum(S(:,run==1),2)/lrun;
fire(:,3) = sum(S(:,run==0),2)/lstand;

firetot(:,1) = sum(S,2);
firetot(:,2) = sum(S(:,run==1),2);
firetot(:,3) = sum(S(:,run==0),2);
%% sum, mean and std of all cells events

meanfire(1,:) = [length(find(fire(:,1))) length(find(fire(:,2))) length(find(fire(:,3)))];
meanfire(2,:) = [mean(fire(:,1)) mean(fire(:,2)) mean(fire(:,3))];
meanfire(3,:) = [std(fire(:,1)) std(fire(:,2)) std(fire(:,3))]/sqrt(size(fire,1));
meanfire(isnan(meanfire)) = 0;
%% Amplitude Comparison

win = 10; % Window size in frames to look for maximum
for i = 1:size(S,1)
    a = C(i,S(i,:) ==1);
    b = find(S(i,:));
    c = 1:size(b,2); 
    d = run(b);
    for j = 1:length(b)  
        if b(j)+win<size(C,2)
            a(j) = max(C(i,b(j):b(j)+win));
            b(j) = b(j)+find(C(i,b(j):b(j)+win)==max(C(i,b(j):b(j)+win)),1)-1;
        end
        if j > 1 && b(j) == b(j-1)
            c(j) = 0;
        end
    end
    a = a(c(c~=0))/Df(i);
    b = b(c(c~=0));
    d = d(c(c~=0));
    if find(a)
        amplitude(i,1) = mean(a);
        maxamplitude(i,1) = max(a);
        if find(d==1)
            amplitude(i,2) = mean(a(d==1)); 
            maxamplitude(i,2) = max(a(d==1));
        else
            amplitude(i,2) = NaN; 
            maxamplitude(i,2) = NaN;
        end
        if find(d==0)
            maxamplitude(i,3) = max(a(d==0));    
            amplitude(i,3) = mean(a(d==0));  
        else
            amplitude(i,3) = NaN; 
            maxamplitude(i,3) = NaN;
        end
    else
        amplitude(i,:) = NaN;
        maxamplitude(i,:) = NaN;
    end
end

for i = 1:3

    meanamplitude(1,i) = mean(amplitude(~isnan(amplitude(:,i)),i));
    meanamplitude(2,i) = std(amplitude(~isnan(amplitude(:,i)),i));
    if find(~isnan(amplitude(:,i)))
        meanamplitude(3,i) = max(maxamplitude(~isnan(maxamplitude(:,i)),i));
    else
        meanamplitude(3,i) = NaN;
    end
end

%% statistical properties

run = find(scn.running==1);
norun = find(scn.running==0);
if length(run)<length(norun)
    norun = norun(1:length(run));
else
    run = run(1:length(norun));
end

C = caim.C./caim.Df;
BsStat = [sum(var(C(:,norun),0,2)) sum(var(C(:,run),0,2))  sum(skewness(C(:,norun)')') sum(skewness(C(:,run)')')  sum(kurtosis(C(:,norun)')') sum(kurtosis(C(:,run)')')];
% x = mean(C,1);
x = sum(C,1);
SumStat = [var(x(norun),0) var(x(:,run),0)  skewness(x(norun),1) skewness(x(:,run))  kurtosis(x(norun),1) kurtosis(x(run),1)];

%%
x = 0:.004:2;
xx = 0:.2:100;
xxx = 0:.001:.5;
a = histcounts(C(:,norun),x,'normalization','probability');
b = histcounts(C(:,run),x,'normalization','probability');
c = histcounts(sum(C(:,norun),1),xx,'normalization','probability');
d = histcounts(sum(C(:,run),1),xx,'normalization','probability');
e = histcounts(mean(C(:,norun),1),xxx,'normalization','probability');
f = histcounts(mean(C(:,run),1),xxx,'normalization','probability');
FireHist = [a;b;c;d;e;f];

%%
if ~isempty(run)
    a = C(:,norun);a = a(:);
    b = C(:,run);b = b(:);
    [~,p1] = kstest2(a,b);
    % ystat = [a(:); b(:)];
    % ygroup(1:length(a(:))) = {'rest'};
    % ygroup(end+1:end+length(b(:))) = {'run'};
    % [p1,~,~] = kruskalwallis(ystat,ygroup,'off');
    a = sum(C(:,norun),1);
    b = sum(C(:,run),1);
    [~,p2] = kstest2(a,b);
    % ystat = [a(:); b(:)];
    % ygroup = {};
    % ygroup(1:length(a(:))) = {'rest'};
    % ygroup(end+1:end+length(b(:))) = {'run'};
    % [p2,~,~] = kruskalwallis(ystat,ygroup,'off');
    a = mean(C(:,norun),1);
    b = mean(C(:,run),1);
    [~,p3] = kstest2(a,b);
    % ystat = [a(:); b(:)];
    % ygroup = {};
    % ygroup(1:length(a(:))) = {'rest'};
    % ygroup(end+1:end+length(b(:))) = {'run'};
    % [p3,~,~] = kruskalwallis(ystat,ygroup,'off');
else
    p1 = NaN; p2 = NaN; p3 = NaN;
end
%% output

fireprop.fire = fire;
fireprop.firetot = firetot;
fireprop.meanfire = meanfire;
fireprop.amplitude = amplitude;
fireprop.meanamplitude = meanamplitude;
fireprop.BsStat = BsStat;
fireprop.SumStat = SumStat;
fireprop.FireHist = FireHist;
fireprop.FireHistp = [p1 p2 p3];
%% plot bargraph with frequencies in different stages

% % 
% bar(meanfire(1,:),'w')
% hold on
% errorbar(meanfire(1,:),meanfire(2,:),'.','linewidth',2,'color',[.6 .6 .6])
% set(gcf,'color',[0 0 0])
% set(gca,'fontsize',20,...
%             'YColor',[1 1 1],...
%             'XColor',[1 1 1],...
%             'Color',[0 0 0],...
%             'XTickLabel',{'Total','Running','Resting'})
% hold off
% 
% % exportfig('Frequency per state')
%%

% histogram(fire(:,2),'BinWidth',1);%,'normalization','probability')
% hold on
% histogram(fire(:,1),'BinWidth',1);%,'normalization','probability')
% hold off
% xtil = string(nbins);
% xtil(end) = string([' >' num2str(nbins(end-1))]); 
% set(gca,'Xlim',[1 30]);%,...
%             'fontsize',20,...
%             'YColor',[1 1 1],...
%             'XColor',[1 1 1],...
%             'Color',[0 0 0],...
%             'XTickLabel',xtil)
% set(gcf,'color',[0 0 0])
% 
% % exportfig('act distribution')

end



