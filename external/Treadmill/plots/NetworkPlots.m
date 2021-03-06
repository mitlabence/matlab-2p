% % network frequency and occurence plots

Cues = load('Z:\Martin Pofahl\Cues.mat');
Baseline = load('Z:\Martin Pofahl\Baseline.mat');
Airpuff = load('Z:\Martin Pofahl\Airpuff.mat');
load('Z:\Martin Pofahl\BigFatCluster.mat');
load('cclustID.mat')

experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};
mouse = {'M103' 'M155' 'M158' 'M194' 'M195' 'M224' 'M226' 'M227' 'M229' 'M234'}; 
plotcol = [0 .3 0;
        0 .6 0;
        0 .9 0;
        .3 0 0;
        .6 0 0;
        .9 0 0;
        0 0 .3;
        0 0 .6;
        0 0 .9];

trial = 3;
ftsz = 8;
lnwd = 1;


%% net freq
figure

x = [1 2 4 5 7 8 10 11 13 14 16 17];
mouseID = Baseline.Fire.mouseID;
yin = Baseline.Network.netfreq(:,2:3);
y = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
yerr = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,2)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,2))==i,:),1);
    yerr(i,:) = nanstd(yin(cell2mat(mouseID(:,2))==i,:),1);
end

y = y(3:5,:);
y = reshape(y',size(y,1)*size(y,2),1);
yerr = yerr(3:5,:);
yerr = reshape(yerr',size(yerr,1)*size(yerr,2),1);

mouseID = Cues.Fire.mouseID;
yin =  Cues.Network.netfreq(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
y1err = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,2)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,2))==i,:),1);
    y1err(i,:) = nanstd(yin(cell2mat(mouseID(:,2))==i,:),1);
end
y1 = reshape(y1',size(y1,1)*size(y1,2),1);
y1err = reshape(y1err',size(y1err,1)*size(y1err,2),1);

y = [y; y1];
yerr = [yerr; y1err];
yerr = yerr./sqrt(max(cell2mat(mouseID(:,4))));

% reps = size(y,1);
% [~,~,stats] = anova2(ystat,reps,'off');
% c = multcompare(stats,'Display','off');
% c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

b = bar(x,y);
hold on

errorbar(x,y,yerr,'.',...
            'Color',[0 0 0],...
            'LineWidth',lnwd)

b.FaceColor = 'flat';
b.LineStyle = 'none';
for i = 1:length(y)/2
    b.CData(2*i-1,:) = [.1 1 .1];
    b.CData(2*i,:) = [0 .5 0];
end

ax = gca;
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XTick = 1.5 :3: 16.5;
ax.XTickLabel = {'Base1' 'Base2' 'Base3' 'Cues1' 'Cues2' 'Cues3'};

title('Frequency of networks')
%% net occurence
figure

x = [1 2 4 5 7 8 10 11 13 14 16 17];
mouseID = Baseline.Fire.mouseID;
yin = Baseline.Network.netoccur(:,2:3);
y = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
yerr = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,2)))
    y(i,:) = nanmean(yin(cell2mat(mouseID(:,2))==i,:),1);
    yerr(i,:) = nanstd(yin(cell2mat(mouseID(:,2))==i,:),1);
end

y = y(3:5,:);
y = reshape(y',size(y,1)*size(y,2),1);
yerr = yerr(3:5,:);
yerr = reshape(yerr',size(yerr,1)*size(yerr,2),1);

mouseID = Cues.Fire.mouseID;
yin =  Cues.Network.netoccur(:,2:3);
y1 = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
y1err = zeros(max(cell2mat(mouseID(:,2))),size(yin,2));
for i = 1:max(cell2mat(mouseID(:,2)))
    y1(i,:) = mean(yin(cell2mat(mouseID(:,2))==i,:),1);
    y1err(i,:) = nanstd(yin(cell2mat(mouseID(:,2))==i,:),1);
end
y1 = reshape(y1',size(y1,1)*size(y1,2),1);
y1err = reshape(y1err',size(y1err,1)*size(y1err,2),1);

y = [y; y1];
yerr = [yerr; y1err];
yerr = yerr./sqrt(max(cell2mat(mouseID(:,4))));

% reps = size(y,1);
% [~,~,stats] = anova2(ystat,reps,'off');
% c = multcompare(stats,'Display','off');
% c(2:4,:) = multcompare(stats,'Estimate','row','Display','off');

b = bar(x,y);
hold on

errorbar(x,y,yerr,'.',...
            'Color',[0 0 0],...
            'LineWidth',lnwd)

b.FaceColor = 'flat';
b.LineStyle = 'none';
for i = 1:length(y)/2
    b.CData(2*i-1,:) = [.1 1 .1];
    b.CData(2*i,:) = [0 .5 0];
end

ax = gca;
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
ax.XTick = 1.5 :3: 16.5;
ax.XTickLabel = {'Base1' 'Base2' 'Base3' 'Cues1' 'Cues2' 'Cues3'};

title('Occurence of networks')
%% Network size cumulative
figure

a1 = [];
a2 = [];
a3 = [];
a4 = [];
a5 = [];
a6 = [];
% a7 = [];
% a8 = [];
% a9 = [];

for i = 1:size(CAIM,2)   
    if ~isempty(CAIM(3,i).A)
        a1 = [a1 CAIM(3,i).network.netnum];
    end
    if ~isempty(CAIM(4,i).A)
        a2 = [a2 CAIM(4,i).network.netnum];
    end
    if ~isempty(CAIM(5,i).A)
        a3 = [a3 CAIM(5,i).network.netnum];
    end
    if ~isempty(CAIM(6,i).A)
        a4 = [a4 CAIM(6,i).network.netnum];
    end
    if ~isempty(CAIM(7,i).A)
        a5 = [a5 CAIM(7,i).network.netnum];
    end
    if ~isempty(CAIM(8,i).A)
        a6 = [a6 CAIM(8,i).network.netnum];
    end
%     if ~isempty(CAIM(9,i).A)
%         a7 = [a7 CAIM(9,i).network.netnum];
%     end
%     if ~isempty(CAIM(10,i).A)
%         a8 = [a8 CAIM(10,i).network.netnum];
%     end
%     if ~isempty(CAIM(11,i).A)
%         a9 = [a9 CAIM(11,i).network.netnum];
%     end
end

ystat = [a1 a2 a3 a4 a5 a6];
ygroup = [];
ygroup(1:length(a1)) = 1;
a = length(a1);
ygroup(a+1:a+length(a2)) = 2;
a = a+length(a2);
ygroup(a+1:a+length(a3)) = 3;
a = a+length(a3);
ygroup(a+1:a+length(a4)) = 4;
a = a+length(a4);
ygroup(a+1:a+length(a5)) = 5;
a = a+length(a5);
ygroup(a+1:a+length(a6)) = 6;
[~,~,stats] = kruskalwallis(ystat,ygroup,'off');
c = multcompare(stats,'Display','off');

[aa1,bb1]  = histcounts(a1, 'Normalization', 'probability'); 
[aa2,bb2]  = histcounts(a2, 'Normalization', 'probability');
[aa3,bb3]  = histcounts(a3, 'Normalization', 'probability');
[aa4,bb4]  = histcounts(a4, 'Normalization', 'probability'); 
[aa5,bb5]  = histcounts(a5, 'Normalization', 'probability');
[aa6,bb6]  = histcounts(a6, 'Normalization', 'probability');
% [aa7,bb7]  = histcounts(a7, 'Normalization', 'probability'); 
% [aa8,bb8]  = histcounts(a8, 'Normalization', 'probability');
% [aa9,bb9]  = histcounts(a9, 'Normalization', 'probability');


plot(bb1(1:end-1)+.5,cumsum(aa1),...
    'linewidth',lnwd,...
    'color',plotcol(1,:))
hold on
plot(bb2(1:end-1)+.5,cumsum(aa2),...
    'linewidth',lnwd,...
    'color',plotcol(2,:))
plot(bb3(1:end-1)+.5,cumsum(aa3),...
    'linewidth',lnwd,...
    'color',plotcol(3,:))
plot(bb4(1:end-1)+.5,cumsum(aa4),...
    'linewidth',lnwd,...
    'color',plotcol(4,:))
plot(bb5(1:end-1)+.5,cumsum(aa5),...
    'linewidth',lnwd,...
    'color',plotcol(5,:))
plot(bb6(1:end-1)+.5,cumsum(aa6),...
    'linewidth',lnwd,...
    'color',plotcol(6,:))
% plot(bb7(1:end-1)+.5,cumsum(aa7),...
%     'linewidth',lnwd,...
%     'color',plotcol(7,:))
% plot(bb8(1:end-1)+.5,cumsum(aa8),...
%     'linewidth',lnwd,...
%     'color',plotcol(8,:))
% plot(bb9(1:end-1)+.5,cumsum(aa9),...
%     'linewidth',lnwd,...
%     'color',plotcol(9,:))

legend(experiment(3:8),'Location','best')

title('Networks size')

%% # of networks over time
figure
hold on

for i = 3:8%1:size(CAIM,1)
    %%
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    netpos = [];
    for j = 1:length(fullCAIM)
        k = fullCAIM(j);      
        netpos = [netpos; CAIM(i,k).network.netpos/1000/60];
    end
    
    [a,b] = histcounts(netpos,120);
    plot(b(1:end-1),cumsum(a),...
        'color',plotcol(i-2,:))
end
legend(experiment(3:8),'Location','best')
title('Networks over time')

%% Inter Event Interval
figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]],...
        'renderer','painters',...
        'visible','on')
    
intEvInt = zeros(6,9,200);
for i = 3:8%1:size(CAIM,1)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
%     intEvInt = [];

    for j = 1:length(fullCAIM)
        k = fullCAIM(j);    
        intEvInt(i,k,:) = CAIM(i,k).network.intEvInt;
        subplot(7,10,(i-3)*10+k)        
%         plot(.1:.1:20,cumsum(CAIM(i,k).network.intEvInt))
        bar(.1:.1:20,CAIM(i,k).network.intEvInt) 
        if i == 3
            title(mouse(k))
        end
        if j == 1
            ylabel([experiment(i) ', # net ev'])
        end
%         if i == 8
%             xlabel('time/s')
%         end
        grid on
        axis tight
    end
    
end

for j = 1:size(CAIM,2)    
    subplot(7,10,6*10+j)
    a = permute(sum(intEvInt(:,j,:),1),[3 1 2]);
%     plot(.1:.1:20,cumsum(a))    
    bar(.1:.1:20,a) 
    if j == 1
        ylabel(['Sum, # net ev'])
    end
    
    xlabel('time/s')
    
    grid on
    axis tight
end

for i = 3:8  
    subplot(7,10,(i-3)*10+10)
    a = permute(sum(intEvInt(i,:,:),2),[3 2 1]);
%     plot(.1:.1:20,cumsum(a))   
    bar(.1:.1:20,a) 
    if i == 3
        title('sum')
    end
    grid on
    axis tight
end
%% place and net activity 

placefld = [];
netplace = [];
plcvctang = [];
for i = 5
    %%
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];    

    for j = 1:length(fullCAIM)
        k = fullCAIM(j);              
        cclust = CAIM(i,k).cclust;
        isplace = cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05;
        plcvctang = [plcvctang ;cclust(isplace,cclustID.plcvctang)];        
        placefld = [placefld; CAIM(i,k).plcfield(isplace,1:150)];
        netplace = [netplace; CAIM(i,k).network.netplace(isplace,:)];
    end  
end
for i = 1:size(netplace,1)
    netplace(i,:) = 1*(smooth(netplace(i,:),1));
end
[plccenter,aa] = sort(plcvctang);
c = find(abs(plccenter)==min(abs(plccenter)));
aa = aa([c:end 1:c-1],:);

% plot(1:40,cumsum(sum(clustsize)/sum(sum(clustsize))),...
%     'LineWidth',lnwd,...
%     'color',plotcol(i-4,:));

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[20 5 [ 2*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 2*8.9])


subplot(1,3,1)
imagesc(placefld(aa,:))
title('Place Act')
subplot(1,3,2)
imagesc(netplace(aa,:))
title('Networks on Place')
subplot(1,3,3)
A = placefld(aa,:).*netplace(aa,:);
A(A>0) = .7;
A(1) = 1;
imagesc(A)
title('1 x 2')

myColorMap = jet(256);
% myColorMap(1,:) = 1;
colormap(myColorMap);
% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\placenet'); 

% hold off

%% succesive place and net activity 
% close all
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 1 [ 2*sqrt(2)*8.9 3*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 3*8.9])

Spacedistbefore = [];
Timedistbefore = [];
Spacedistafter = [];
Timedistafter = [];
Spaceabsafter = [];
Spaceabsbefore = [];

for i = 8
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];    

    for j = 1:length(fullCAIM)
        k = fullCAIM(j);   
        spacedistbefore = [];
        timedistbefore = [];
        spaceabsafter = [];
        spacedistafter = [];
        timedistafter = [];
        spaceabsbefore = [];
        
        if isfield(CAIM(i,k).network,'fieldnet') && ~isempty(CAIM(i,k).network.fieldnet)        
            fieldnet = CAIM(i,k).network.fieldnet;
            
            for ii = 1:size(fieldnet,3)

                % Distance before place field crossing
                fntemp = fieldnet(:,:,ii);
                fntemp = fntemp(~isnan(fntemp(:,5)),:);
                isdouble = false(size(fntemp,1),1);
                for jj = 2:size(fntemp,1)
                    isdouble(jj) = fntemp(jj,3)==fntemp(jj-1,3);
                end
                fntemp = fntemp(~isdouble,:);
                spacedistbefore = [spacedistbefore; abs(fntemp(:,6)-fntemp(:,2))];
                spaceabsbefore = [spaceabsbefore; fntemp(:,7) fntemp(:,8)];
                timedistbefore = [timedistbefore; abs(fntemp(:,5)-fntemp(:,1))]; 

                % Distance after place field crossing
                fntemp = fieldnet(:,:,ii);
                fntemp = fntemp(~isnan(fntemp(:,3)),:);
                isdouble = false(size(fntemp,1),1);
                for jj = 1:size(fntemp,1)-1
                    isdouble(jj) = fntemp(jj,5)==fntemp(jj+1,5);
                end
                fntemp = fntemp(~isdouble,:);
                spacedistafter = [spacedistafter; abs(fntemp(:,4)-fntemp(:,2))]; 
                spaceabsafter = [spaceabsafter; [fntemp(:,7) fntemp(:,8)]]; 
                timedistafter = [timedistafter; abs(fntemp(:,3)-fntemp(:,1))]; 
            end
            inround = spacedistbefore<1500;
            spaceabsbefore = spaceabsbefore(inround,:);
            timedistbefore = timedistbefore(inround);
            spacedistbefore = spacedistbefore(inround);
            inround = spacedistafter<1500;
            spaceabsafter = spaceabsafter(inround,:);
            timedistafter = timedistafter(inround);
            spacedistafter = spacedistafter(inround);
            
            
            Spacedistbefore = [Spacedistbefore; spacedistbefore];           
            Timedistbefore = [Timedistbefore; timedistbefore];
            Spaceabsbefore = [Spaceabsbefore; spaceabsbefore];
            Spacedistafter = [Spacedistafter; spacedistafter];
            Timedistafter = [Timedistafter; timedistafter];
            Spaceabsafter = [Spaceabsafter; spaceabsafter];
            
        end

        subplot(length(fullCAIM)+1,4,1+(j-1)*4)
        histogram(spacedistbefore/10,0:10:150,'FaceColor',[0 0 1])
        ylabel(mouse(k))
        if j == 1;title('Distance field to last net');end
        subplot(length(fullCAIM)+1,4,2+(j-1)*4)
        histogram(timedistbefore/1000,0:10:250,'FaceColor',[1 0 0])
        if j == 1;title('Time field to last net');end
        subplot(length(fullCAIM)+1,4,3+(j-1)*4)
        histogram(spacedistafter/10,0:10:150,'FaceColor',[0 0 1])
        if j == 1;title('Distance field to next net');end
        subplot(length(fullCAIM)+1,4,4+(j-1)*4)
        histogram(timedistafter/1000,0:10:250,'FaceColor',[1 0 0])
        if j == 1;title('Time field to next net');end
    end  
end

subplot(length(fullCAIM)+1,4,(length(fullCAIM))*4+1)
histogram(Spacedistbefore/10,0:10:150,'FaceColor',[0 0 1],'LineWidth',2)
% histogram(Spaceabsbefore/10,0:10:150,'FaceColor',[0 0 1],'LineWidth',2)
xlabel('distance/cm')
ylabel('sum')

subplot(length(fullCAIM)+1,4,(length(fullCAIM))*4+2)
histogram(Timedistbefore/1000,0:10:250,'FaceColor',[1 0 0],'LineWidth',2)
xlabel('delay/s')
ylim([0 10])
text(100,8,['sum = ' num2str(length(Timedistbefore))])

subplot(length(fullCAIM)+1,4,(length(fullCAIM))*4+3)
histogram(Spacedistafter/10,0:10:150,'FaceColor',[0 0 1],'LineWidth',2)
% histogram(Spaceabsafter/10,0:10:150,'FaceColor',[0 0 1],'LineWidth',2)
xlabel('distance/cm')

subplot(length(fullCAIM)+1,4,(length(fullCAIM))*4+4)
histogram(Timedistafter/1000,0:10:250,'FaceColor',[1 0 0],'LineWidth',2)
xlabel('delay/s')
ylim([0 10])
text(100,8,['sum = ' num2str(length(Timedistafter))])

% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\field-net-delay BASELINE'); 

%% scatter plot of Place field and Network position

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 1 [ sqrt(2)*8.9 8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [sqrt(2)*8.9 8.9])

equality = 200;

subplot(1,2,1)
a = Spaceabsbefore;
scatter(a(:,2)/10,a(:,1)/10)
hold on
scatter(a(:,2)/10+150,a(:,1)/10+150)
scatter(a(:,2)/10,a(:,1)/10+150)
scatter(a(:,2)/10+150,a(:,1)/10)
plot([0 equality],[0 equality],'--')
xlim([0 equality])
ylim([0 equality])
xlabel('Network act location /cm')
ylabel('Place field center /cm')
title('Network before field')

subplot(1,2,2)
a = Spaceabsafter;
scatter(a(:,2)/10,a(:,1)/10)
hold on
scatter(a(:,2)/10+150,a(:,1)/10+150)
scatter(a(:,2)/10,a(:,1)/10+150)
scatter(a(:,2)/10+150,a(:,1)/10)
plot([0 equality],[0 equality],'--')
xlim([0 equality])
ylim([0 equality])
xlabel('Network act location /cm')
ylabel('Place field center /cm')
title('Field before network')

% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\field-net scatter');

%% polar plot of Place field and network activity position
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 1 [ sqrt(2)*8.9 8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [sqrt(2)*8.9 8.9])

numbin = 150;
Alpha = linspace(0,2*pi,numbin);
rho = ones(1,length(Alpha));%(1:length(theta))/length(theta);

a = Spaceabsbefore;
aa = zeros(size(a,1),1);
d = zeros(size(a,1),1);
for i = 1:size(a,1)   % loop over cells                    
    aa(i,3) = round((a(i,1)-a(i,2))/10);
    if aa(i,3)<0;aa(i,3) = aa(i,3)+150;end
    if aa(i,3)==0;aa(i,3) = 1;end
    aa(i,:) =   [Alpha(aa(i,3)),rho(aa(i,3)), aa(i,3)];
    d(i) = exp(Alpha(aa(i,3))*1j);      
 
end
prefvec = sum(d)/length(d); % place coding vector

subplot(1,2,1)
h = polarplot(aa(:,1),aa(:,2),'o');
hh = get(h,'color');
hold on
polarplot(prefvec,'*','Color',hh);
polarplot([0 real(-1j*log(prefvec))],[0 abs(prefvec)],'Color',hh);
ax = gca;
% ax.ThetaTickLabel = {};
ax.RGrid = 'off';
ax.RTickLabel = {};
hold off
title('Network before field')

a = Spaceabsafter;
aa = zeros(size(a,1),1);
d = zeros(size(a,1),1);
for i = 1:size(a,1)   % loop over cells                    
    aa(i,3) = round((a(i,1)-a(i,2))/10);
    if aa(i,3)<0;aa(i,3) = aa(i,3)+150;end
    if aa(i,3)==0;aa(i,3) = 1;end
    aa(i,:) =   [Alpha(aa(i,3)),rho(aa(i,3)), aa(i,3)];
%     aa(i,:) =   [Alpha(aa(i,3)),i, aa(i,3)];
    
    d(i) = exp(Alpha(aa(i,3))*1j);      
 
end
prefvec = sum(d)/length(d); % place coding vector

subplot(1,2,2)
h = polarplot(aa(:,1),aa(:,2),'o');
hh = get(h,'color');
hold on
polarplot(prefvec,'*','Color',hh);
polarplot([0 real(-1j*log(prefvec))],[0 abs(prefvec)],'Color',hh);
ax = gca;
% ax.ThetaTickLabel = {};
ax.RGrid = 'off';
ax.RTickLabel = {};
hold off
title('Field before network')
% print(gcf, '-dpdf', 'field-net polar');
% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\field-net polar');

%% network place with respect to place field center for all network events
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

numbin = 150;
Alpha = linspace(0,2*pi,numbin)';
rho = ones(length(Alpha),1);%(1:length(theta))/length(theta);
% Netang = [];

for i = 3:8
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];    
    Netang = [];
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        if isfield(CAIM(i,k).network,'cnet') && ~isempty(CAIM(i,k).network.cnet)            
            cnet = CAIM(i,k).network.cnet;
            netang = [];
            for ii = 1:size(cnet,1)
                netspot = cnet{ii,5};
                netspot = ceil(netspot/10);
                netspot = netspot-ceil(cnet{ii,2}/10);
                netspot(netspot<=0) = netspot(netspot<=0)+150;
                netang =  [netang; Alpha(netspot),rho(netspot), netspot];
            end
            Netang = [Netang; netang];
            if ~isempty(netang)
                subplot(6,10,(i-3)*10+k) 
                d = exp(Alpha(netang(:,3))*1j); 
                prefvec = sum(d)/length(d);  
                
                h = polarplot(netang(:,1),netang(:,2),'o');
                hh = get(h,'color');
                hold on
                polarplot(prefvec,'*','Color',hh);
                polarplot([0 real(-1j*log(prefvec))],[0 abs(prefvec)],'Color',hh);
                ax = gca;
                ax.ThetaTickLabel = {};
                ax.RGrid = 'off';
                ax.RTickLabel = {};
                hold off
                if i == 3;title(mouse(j));end
%                 if j == 1;ylabel(experiment(i));end
            end
        end
    end
    subplot(6,10,(i-3)*10+10) 
    d = exp(Alpha(Netang(:,3))*1j); 
    prefvec = sum(d)/length(d);  
    h = polarplot(Netang(:,1),Netang(:,2),'o');
    hh = get(h,'color');
    hold on
    polarplot(prefvec,'*','Color',hh);
    polarplot([0 real(-1j*log(prefvec))],[0 abs(prefvec)],'Color',hh);
    ax = gca;
    ax.ThetaTickLabel = {};
    ax.RGrid = 'off';
    ax.RTickLabel = {};
    hold off
    if i == 3;title('mean');end
end

% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\field-net ind polar');

%%

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 2*8.9])

d = exp(Alpha(Netang(:,3))*1j); 
prefvec = sum(d)/length(d);  
h = polarplot(Netang(:,1),Netang(:,2),'o');

hh = get(h,'color');
hold on
polarplot(prefvec,'*','Color',hh);
polarplot([0 real(-1j*log(prefvec))],[0 abs(prefvec)],'Color',hh);
ax = gca;
% ax.ThetaTickLabel = {};
ax.RGrid = 'off';
ax.RTickLabel = {};
hold off

%% participants plot

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 2*8.9])

num = [3:5 8];
% bin = [0:.01:1];
bin = [0:1:100];
NetCum = NaN(size(CAIM,2),length(bin)-1,max(num));
NetPlCum = NaN(size(CAIM,2),length(bin)-1,max(num));    
NetSpCum = NaN(size(CAIM,2),length(bin)-1,max(num));
NetCumSum = NaN(length(bin)-1,max(num));
NetPlCumSum = NaN(length(bin)-1,max(num));
NetSpCumSum = NaN(length(bin)-1,max(num));
a = [];
b = [];
c = [];
for i = num%1:size(CAIM,1)   
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    thresh = 0;
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).network) || length(CAIM(i,fullCAIM(j)).network.netpos)<thresh
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    if i>6
        a = [];
        b = [];
        c = [];
    end
    for j = 1:length(fullCAIM)
        
        k = fullCAIM(j);  
        cclusttemp = CAIM(i,k).cclust;
        cclusttemp(:,cclustID.netprob) = cclusttemp(:,cclustID.netprob)*length(CAIM(i,fullCAIM(j)).network.netpos);
        isnet = true(size(cclusttemp,1),1);%
%         isnet = cclusttemp(:,cclustID.netprob)>0;
        isplace = cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05;
        isspeed = cclusttemp(:,cclustID.speedcorr)>.9 & cclusttemp(:,cclustID.speedcorrp2)<=.05;
        CAIM(i,k).network;
        a = [a; cclusttemp(isnet & ~isplace & ~isspeed,cclustID.netprob)];
        b = [b; cclusttemp(isnet & isplace,cclustID.netprob)];
        c = [c; cclusttemp(isnet & isspeed,cclustID.netprob)];
        aa = histcounts(cclusttemp(isnet,cclustID.netprob),bin);
        NetCum(k,:,i) = cumsum(aa)/sum(aa);
        aa = histcounts(cclusttemp(isnet&isplace,cclustID.netprob),bin);
        NetPlCum(k,:,i) = cumsum(aa)/sum(aa);
        aa = histcounts(cclusttemp(isnet&isspeed,cclustID.netprob),bin);
        NetSpCum(k,:,i) = cumsum(aa)/sum(aa);
    end
    aa = histcounts(a,bin);
    NetCumSum(:,i) = cumsum(aa)/sum(aa);
    aa = histcounts(b,bin);
    NetPlCumSum(:,i) = cumsum(aa)/sum(aa);
    aa = histcounts(c,bin);
    NetSpCumSum(:,i) = cumsum(aa)/sum(aa);
end

subplot(2,3,1)
plot(bin(1:end-1),(NetCum(:,:,5)),'color',[.5 .5 .5])
hold on
plot(bin(1:end-1),NetCumSum(:,5),'Linewidth',2)
hold off
title('All cells')
subplot(2,3,2)
plot(bin(1:end-1),(NetPlCum(:,:,5)),'color',[.5 .5 .5])
hold on
plot(bin(1:end-1),NetPlCumSum(:,5),'Linewidth',2)
hold off
title('Place cells')
subplot(2,3,3)
plot(bin(1:end-1),(NetSpCum(:,:,5)),'color',[.5 .5 .5])
hold on
plot(bin(1:end-1),NetSpCumSum(:,5),'Linewidth',2)
hold off
title('Speed cells')
subplot(2,3,4)
plot(bin(1:end-1),(NetCum(:,:,8)),'color',[.5 .5 .5])
hold on
plot(bin(1:end-1),NetCumSum(:,8),'Linewidth',2)
hold off
title('All cells')
subplot(2,3,5)
plot(bin(1:end-1),(NetPlCum(:,:,8)),'color',[.5 .5 .5])
hold on
plot(bin(1:end-1),NetPlCumSum(:,8),'Linewidth',2)
hold off
title('Place cells')
subplot(2,3,6)
plot(bin(1:end-1),(NetSpCum(:,:,8)),'color',[.5 .5 .5])
hold on
plot(bin(1:end-1),NetSpCumSum(:,8),'Linewidth',2)
hold off
title('Speed cells')

% print(gcf, '-dpdf', 'C:\Users\Admin\Dropbox\Dentate in-vivo Project\Net Prob');
%% Network and cell number

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2*sqrt(2)*8.9 2*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2*sqrt(2)*8.9 2*8.9])

num = [1:5];
hold on
for i = num%1:size(CAIM,1)   
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(i,fullCAIM(j)).network) || ~isfield(CAIM(i,fullCAIM(j)).network,'netraster')
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
    
    for j = 1:length(fullCAIM)
        
        k = fullCAIM(j);
            a = size(CAIM(i,k).network.netraster);
            scatter(a(1),a(2))

    end
    
    
    
end