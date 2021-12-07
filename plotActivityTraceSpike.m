function plotActivityTraceSpike( inferred, filtered, raw, spike, plotControl)

% this code plots the activity trace with spike

if isfield(plotControl,'frameRate') frameRate=plotControl.frameRate;
else frameRate=1; end

if isfield(plotControl,'normalization') normalization=plotControl.normalization;
else normalization=1; end

if isfield(plotControl,'sep') sep=plotControl.sep;
else sep=1.2; end

if isfield(plotControl,'displayLabel') displayLabel=plotControl.displayLabel;
else displayLabel=1; end

if isfield(plotControl,'plotInferred') plotInferred=plotControl.plotInferred;
else plotInferred=1; end

if isfield(plotControl,'plotFiltered') plotFiltered=plotControl.plotFiltered;
else plotFiltered=1; end

if isfield(plotControl,'plotRaw') plotRaw=plotControl.plotRaw;
else plotRaw=0; end

if isfield(plotControl,'plotSpike') plotSpike=plotControl.plotSpike;
else plotSpike=1; end

if isfield(plotControl,'rollingView') rollingView=plotControl.rollingView;
else rollingView=0; end

t=1/frameRate*(1:size(spike,2));
index=find(spike==0);
spike(index)=NaN;
K=size(inferred,1);

figure;
hold on;

if normalization==1
    for idx=1:K
        if plotRaw  plot(t,raw(idx,:)/max(raw(idx,:))-idx*sep,'color',[0.8 0.8 0.8],'linewidth',2); end
        if plotFiltered  plot(t,filtered(idx,:)/max(filtered(idx,:))-idx*sep,'color',[1 0.7 0.7],'linewidth',1.5); end
        if plotInferred  plot(t,inferred(idx,:)/max(inferred(idx,:))-idx*sep,'color',[0 0 1],'linewidth',1); end
    end
    ylabel('Normalized \DeltaF/F');
else
    for idx=1:K
        if plotRaw  plot(t,raw(idx,:)-idx*sep,'color',[0.8 0.8 0.8],'linewidth',2); end
        if plotFiltered  plot(t,filtered(idx,:)-idx*sep,'color',[1 0.7 0.7],'linewidth',1.5); end
        if plotInferred  plot(t,inferred(idx,:)-idx*sep,'color',[0 0 1],'linewidth',1); end
    end    
    ylabel('\DeltaF/F');
end
xlabel('Time (s)');

if plotSpike
    for idx=1:K
        scatter(t,spike(idx,:)-1-idx*sep-0.2, '.', 'linewidth',2,'MarkerEdgeColor',[0 0 0]);
    end
end

if displayLabel
    for idx=1:K
        text(max(t)*1.02,-idx*sep, num2str(idx));
    end
end

h=gca;
if ~rollingView 
    axis(h,[0 max(t) -sep*(K+1) sep]);
else
    idx=1;
    displayNum=20;
    while idx+displayNum<K
        axis(h,[0 max(t) -sep*(idx+displayNum) sep*(-idx+1)]);
        idx=idx+displayNum;
        input('Press ENTER to continue...');
    end
    axis(h,[0 max(t) -sep*(idx+displayNum) sep*(-idx+1)]);
    input('Press ENTER to continue...');
    axis(h,[0 max(t) -sep*(K+1) sep]);    
end

