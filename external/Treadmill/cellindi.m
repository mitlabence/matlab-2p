function cellindi(CAIMcorr,mouse,expID)

experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};
%%
load('cclustID.mat','cclustID');
meanf = cclustID.meanf;             % mean frequency of cell firing

plcfld =  cclustID.plcfld;          % placefield center
% plcfldp = cclustID.plcfldp;       % placefield p value
plclength = cclustID.plclength;     % size of place field 
plcvang = cclustID.plcvctang;       % place coding vector angle
plcvct = cclustID.plcvct;           % place coding vector length
plcvctp = cclustID.plcvctp;         % place coding vector p value
runperc = cclustID.runperc;         % percentage of running associated cells
runnerp = cclustID.runnerp;         % p value running

nstim = cclustID.nstim;             % # responded stimuli
airpp = cclustID.airpp;             % p value number of response events
stimang = cclustID.stimang;

mpp = cclustID.mpp;                 % mean MPP input

netperc = cclustID.netperc;         % percentage of network related firing
netprob = cclustID.netprob;         % netprob
corrstimposmean = cclustID.corrstimposmean;           % mean of sig pos place correlations
corrstimnegmean = cclustID.corrstimnegmean;     % 36 mean of sig neg place correlations
corrplaceposmean = cclustID.corrplaceposmean;     % 24 mean of sig pos place correlations
corrposmean = cclustID.corrposmean;     % 18 mean of sig pos correlations

followcell = cclustID.followcell;   % logical if cell was in FOV of all experiments
%%
emptyCAIM = [];
for i = 1:length(expID)
    if isempty(CAIMcorr(expID(i),mouse).A)
        emptyCAIM = [emptyCAIM i];
    end
end
expID(emptyCAIM) = [];
cclust = cell(length(expID),1);
netraster = cell(length(expID),1);
netnum = zeros(length(expID)+1,1);
isplace = cell(length(expID),1);
isstim = cell(length(expID),1);
isboth = cell(length(expID),1);
isthere = cell(length(expID),1);

for i = 1:length(expID)
    trial = expID(i);
    netraster{i} = CAIMcorr(trial,mouse).network.netraster;
    netnum(i+1) = size(netraster{i},2)+netnum(i);
    cclusttemp = CAIMcorr(trial,mouse).cclust;
%     active field sorting
%     netraster{i} = netraster{i}(CAIMcorr(trial,mouse).inField==1,:);
%     cclusttemp = cclusttemp(CAIMcorr(trial,mouse).inField==1,:);
    netraster{i}(CAIMcorr(trial,mouse).inField==0,:) = 0;
    cclusttemp(CAIMcorr(trial,mouse).inField==0,:) = 0;

    isplace{i} = (cclusttemp(:,plcvct)>0 & cclusttemp(:,plcvctp)<=.05) & ~(cclusttemp(:,nstim)>2 & cclusttemp(:,airpp)<=.05);
    isstim{i}  = ~(cclusttemp(:,plcvct)>0 & cclusttemp(:,plcvctp)<=.05) & (cclusttemp(:,nstim)>2 & cclusttemp(:,airpp)<=.05);
    isboth{i}  = (cclusttemp(:,plcvct)>0 & cclusttemp(:,plcvctp)<=.05) & (cclusttemp(:,nstim)>2 & cclusttemp(:,airpp)<=.05);
    isthere{i} = ~isnan(cclusttemp(:,1));
    cclust{i} = cclusttemp;
end

%% fill up individual indicators so that they all have smae length

for i = 1:length(expID-1)
    netraster{i}(size(netraster{end},1),:) = 0;
    isplace{i}(size(netraster{end},1)) = 0;
    isstim{i}(size(netraster{end},1)) = 0;
    isboth{i}(size(netraster{end},1)) = 0;
    isthere{i}(size(netraster{end},1)) = 0;
    cclust{i}(size(netraster{end},1),1) = 0;
end

%%
cclusttemp = cell2mat(cclust(1));
for i = 2:length(cclust)
    cclusttemp = cat(3,cclusttemp,cell2mat(cclust(i)));
end
cclust = cclusttemp;

% transfer the individual cells to one matrix with aditional indicator in
% first coloumn
isplace = cell2mat(isplace');
isplace(:,2:end+1) = isplace; 
isplace(:,1) = sum(isplace(:,2:end),2);
isplace(isplace>1) = 1;

isstim = cell2mat(isstim');
isstim(:,2:end+1) = isstim; 
isstim(:,1) = sum(isstim(:,2:end),2);
isstim(isstim>1) = 1;

isboth = cell2mat(isboth');
isboth(:,2:end+1) = isboth; 
isboth(:,1) = sum(isboth(:,2:end),2);
isboth(isboth>1) = 1;


isboth(isstim(:,1)==1 & isplace(:,1)==1,1) = 1;
isstim(isboth(:,1)==1,1)=0;
isstim(isplace(:,1)==1,1) = 0;
isplace(isboth(:,1)==1,1)=0;
isplace(isstim(:,1)==1,1)=0;


isthere = cell2mat(isthere');
isthere(:,2:end+1) = isthere;

%% read out individual features of our friends

cellID = find(isplace(:,1));
cclusttemp = cclust(cellID,:,:);
cclusttemp = permute(cclusttemp,[1 3 2]);


for i = 1:size(cclusttemp,1)
    %%
    h = figure('color',[1 1 1],...
        'position',[500 50 1.5*[594 420]],...
        'renderer','painters',...
        'visible','on');
    %%
    for j = 1:size(cclusttemp,2)
        %%
        cm = CAIMcorr(expID(j),mouse).cm;
        if cellID(i)<=size(cm,1)
            hh=subplot(6,length(expID),j);
            [d1,d2] = size(CAIMcorr(expID(j),mouse).Cn);
            sx = 10;%min([floor(d1/2),floor(d2/2)]);
            int_x = round(cm(cellID(i),1)) + (-(sx-1):sx);
            if int_x(1)<1
                int_x = int_x + 1 - int_x(1);
            end
            if int_x(end)>d1
                int_x = int_x - (int_x(end)-d1);
            end
            int_y = round(cm(cellID(i),2)) + (-(sx-1):sx);
            if int_y(1)<1
                int_y = int_y + 1 - int_y(1);
            end
            if int_y(end)>d2
                int_y = int_y - (int_y(end)-d2);
            end
            
            B = reshape(CAIMcorr(expID(j),mouse).B(:,cellID(i)),d1,d2);
            B = B(int_x,int_y);
            if find(B,1)
                imagesc(int_x,int_y,B);
                colormap(hh,'gray')
%                 axis off
            else
                axis off
            end
        else
            subplot(6,length(expID),j);
            axis off
        end
        if j == 1
            a = cclusttemp(i,:,cclustID.expID)>0;
            title(['Cell ' num2str(cellID(i)) ', M' num2str(cclusttemp(i,find(a,1),cclustID.expID))]);
        end
        
    end
    
    %%
    
    y = cclusttemp(i,:,plcvang);
    yy = cclusttemp(i,:,cclustID.plcvct);
    yp = cclusttemp(i,:,cclustID.plcvctp);
    for j = 1:length(y)
        subplot(6,length(y),j+length(y))
        polarplot([0 y(j)],[0 yy(j)]) 
        ax = gca;
        ax.ThetaTickLabel = {};
        ax.RLim = [0 1];
        ax.RGrid = 'off';
        ax.RTickLabel = {};
        if j == 1
            if yy(j) > 0 && yp(j)<.05
                title(['plc ang *']);
            else
                title(['plc ang']);
            end
        elseif yy(j) >0 && yp(j)<.05
               title('*')
        end
    end
    
    y = cclusttemp(i,:,cclustID.stimang);
    yy = y~=0;
    yp = cclusttemp(i,:,cclustID.airpp)<.05 & cclusttemp(i,:,cclustID.nstim)>1;
    for j = 1:length(y)
        subplot(6,length(y),j+2*length(y))
        polarplot([0 y(j)],[0 yy(j)]) 
        ax = gca;
        ax.ThetaTickLabel = {};
        ax.RGrid = 'off';
        ax.RTickLabel = {};        
        if j == 1
            if yp(j)
                title(['stim ang *']);
            else
                title(['stim ang']);
            end
        elseif yp(j)
               title('*')
        end
    end 
    
%     subplot(12,1,5)
%     y = cclusttemp(i,:,cclustID.corrpossum);
%     bar(1:length(y),y)
%     ax = gca;
%     ax.XLim = [0.5 length(y)+.5];
%     ax.XTickLabel = [];
%     title('net Prob')
%     
    subplot(10,1,6)
    y = cclusttemp(i,:,cclustID.netprob);
    imagesc(y)
    colormap(jet)
    caxis([0 1])
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    title('Netprob')
    
    subplot(10,1,7)
    y = cclusttemp(i,:,cclustID.corrpossum);
    bar(1:length(y),y)
    ax = gca;
    ax.XLim = [0.5 length(y)+.5];
    ax.XTickLabel = [];
    title('# correlated cells')
    
    subplot(10,1,8)
    y = cclusttemp(i,:,cclustID.corrposmean);
    imagesc(y)
    colormap(jet)
    caxis([0 1])
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    title('mean r')
    
    subplot(10,1,9)
    y = cclusttemp(i,:,cclustID.corrplacepossum);
    bar(1:length(y),y)
    ax = gca;
    ax.XLim = [0.5 length(y)+.5];
    ax.XTickLabel = [];
    title('# correlated PCs')
    
    subplot(10,1,10)
    y = cclusttemp(i,:,cclustID.corrplaceposmean);
    imagesc(y)
    colormap(jet)
    ax = gca;
    ax.XTick = 1:length(y);
    ax.YTick = [];
    ax.XTickLabel = experiment(expID);
    caxis([0 1])
    title('mean r')
    
%     subplot(14,1,11)
%     y = cclusttemp(i,:,cclustID.corrstimpossum);
%     bar(1:length(y),y)
%     ax = gca;
%     ax.XLim = [0.5 length(y)+.5];
%     ax.XTickLabel = [];
%     title('stim coor')
% 
%     subplot(14,1,12)
%     y = cclusttemp(i,:,cclustID.corrstimposmean);
%     imagesc(y)
%     colormap(jet)
%     ax = gca;
%     ax.XTick = [];
%     ax.YTick = [];
%     caxis([0 1])
    %%
%     printpdf(['M' num2str(cclusttemp(i,find(a,1),cclustID.expID)) ', Cell ' num2str(cellID(i)) ])
%     delete(h)
end


end




