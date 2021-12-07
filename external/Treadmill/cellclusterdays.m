function [ClustProp,partnerhist] = cellclusterdays(CAIMcorr,mouse,expID)
% pairwise analysis if cells stay cluster partners on consecutive days to
% calculate over all probability

experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};
load('cclustID.mat','cclustID')
emptyCAIM = [];
for i = 1:length(expID)
    if isempty(CAIMcorr(expID(i),mouse).A)
        emptyCAIM = [emptyCAIM i];
    end
end
expID(emptyCAIM) = [];

%% read out daily clusters and  go through every cluster for pairwise analysis

ClustProp = cell(0,0);
partnerhist = zeros(length(expID)-1,10);
partnerhist2 = zeros(length(expID)-1,10);

for i = 1:length(expID)-1
    %%
    trial = expID(i);
    cclust1 = CAIMcorr(trial,mouse).cclust;
%     cclust1(CAIMcorr(trial,mouse).inField==0,:) = NaN;
    cclust2 = CAIMcorr(expID(i+1),mouse).cclust;
%     cclust2(CAIMcorr(expID(i+1),mouse).inField==0,:) = NaN;
    isplace = (cclust1(:,cclustID.plcvct)>0 & cclust1(:,cclustID.plcvctp)<=.05);
    isstim = (cclust1(:,cclustID.nstim)>2 & cclust1(:,cclustID.airpp)<=.05);
    isplace1 = isplace & ~isstim;
    isstim1  = ~isplace & isstim;
    isboth1 = isplace & isstim;
    isplace = (cclust2(:,cclustID.plcvct)>0 & cclust2(:,cclustID.plcvctp)<=.05);
    isstim = (cclust2(:,cclustID.nstim)>2 & cclust2(:,cclustID.airpp)<=.05);
    isplace2 =  isplace & ~isstim;
    isstim2  = ~isplace & isstim;
    isboth2 = isplace & isstim;
    
    partnerprob = NaN(size(cclust1,1),1);
    stblclust = cell(7,max(cclust1(:,cclustID.clust)));
    
    for j = 1:size(cclust1,1)
        if ~isnan(cclust1(j,cclustID.clust)) && ~isnan(cclust2(j,cclustID.clust))
            % number of partners
            partners = find(cclust1(:,cclustID.clust) == cclust1(j,cclustID.clust));
            partners(partners==j) = [];
            % number of new partners
            newpartners = find(cclust2(:,cclustID.clust) == cclust2(j,cclustID.clust));
            newpartners(newpartners == j) = [];
            samepartner = zeros(size(partners));
            for k = 1:length(partners)
                if find(newpartners == partners(k))
                    samepartner(k) = 1;
                end
            end
            if isempty(partners)
                partnerprob(j) = 0;
            else
                partnerprob(j) = sum(samepartner);%/length(samepartner);
            end
        elseif ~isnan(cclust1(j,cclustID.clust)) && isnan(cclust2(j,cclustID.clust))
            partnerprob(j) = 0;
        end
        if ~isnan(partnerprob(j))%partnerprob(j)>0
            if isempty(find(cclust1(j,cclustID.clust)==cell2mat(stblclust(1,:)),1))
                stblclust{1,cclust1(j,cclustID.clust)} = cclust1(j,cclustID.clust);
                if partnerprob(j)>0 
                    stblclust{2,cclust1(j,cclustID.clust)} = cclust2(j,cclustID.clust);
                end
                stblclust{3,cclust1(j,cclustID.clust)} = length(find(cclust1(j,cclustID.clust)==cclust1(:,cclustID.clust)));
                stblclust{4,cclust1(j,cclustID.clust)} = 1;
            elseif partnerprob(j)>0               
                stblclust{2,cclust1(j,cclustID.clust)} = [stblclust{2,cclust1(j,cclustID.clust)} cclust2(j,cclustID.clust)];
                stblclust{4,cclust1(j,cclustID.clust)} = stblclust{4,cclust1(j,cclustID.clust)} +1;
            end
        end        
    end
    partnerprob(~isnan(partnerprob)) = partnerprob(~isnan(partnerprob))+1;
    partnerhist(i,:) = histcounts(partnerprob,1:1:11);%,'normalization','probability'); %(partnerprob>0)
    
    for j = 1:size(stblclust,2)
        if ~isempty(stblclust{1,j})
            stblclust{5,j} = sum(isplace1(cclust1(:,cclustID.clust) == stblclust{1,j}));
            stblclust{6,j} = sum(isstim1(cclust1(:,cclustID.clust) == stblclust{1,j}));
            stblclust{7,j} = sum(isboth1(cclust1(:,cclustID.clust) == stblclust{1,j}));
        end
    end
    
    %% mean r-value per day and cluster

    netcorr = CAIMcorr(trial,mouse).network.netcorr;
    netcorr(1:1+size(netcorr,1):end) = NaN;
    netcorr2 = CAIMcorr(expID(i+1),mouse).network.netcorr;
    netcorr2(1:1+size(netcorr2,1):end) = NaN;
%     e = 0;
%     A = [];
    for j = 1:size(stblclust,2)
        %%
        k = max(cclust1(:,cclustID.clust))+1-j;
        b = cclust1(:,cclustID.clust) == k;
%         A(e+1:e+sum(b),e+1:e+sum(b)) = netcorr(b,b);
%         e = e+sum(b);
        stblclust{8,k} = nanmean(nanmean(netcorr(b,b)));      
        stblclust{9,k} = nanmean(nanmean(netcorr2(b,b)));
        c = zeros(1,size(stblclust,2)-1);
        for kk = 1:size(stblclust,2)
            if kk~=k
                bb = cclust1(:,cclustID.clust) == kk;
                c(kk) = nanmean(nanmean(netcorr(b,bb))); 
            end
        end
        stblclust{10,k} = c;
    end
%     figure
%     plot(cell2mat(stblclust(8,:)))
%     hold on
%     plot(cell2mat(stblclust(9,:)))
%     hold off
    %%
    partnerprob2 = NaN(size(cclust2,1),1);
    stblclust2 = cell(7,max(cclust2(:,cclustID.clust)));
    cclust1(end+1:size(cclust2,1),:) = NaN;
    isstim1(end+1:size(isstim2,1)) = 0;
    isplace1(end+1:size(isplace2,1)) = 0;
    isboth1(end+1:size(isboth2,1)) = 0;
    
    for j = 1:size(cclust2,1)
        if  ~isnan(cclust2(j,cclustID.clust)) && j <= size(cclust1,1)  && ~isnan(cclust1(j,cclustID.clust)) 
            % number of partners
            partners = find(cclust2(:,cclustID.clust) == cclust2(j,cclustID.clust));
            partners(partners==j) = [];
            % number of new partners
            newpartners = find(cclust1(:,cclustID.clust) == cclust1(j,cclustID.clust));
            newpartners(newpartners == j) = [];
            samepartner = zeros(size(partners));
            for k = 1:length(partners)
                if find(newpartners == partners(k))
                    samepartner(k) = 1;
                end
            end
            if isempty(partners)
                partnerprob2(j) = 0;
            else
                partnerprob2(j) = sum(samepartner)/length(samepartner);
            end
        elseif ~isnan(cclust2(j,cclustID.clust)) && (j > size(cclust1,1) || isnan(cclust1(j,cclustID.clust)))
            partnerprob2(j) = 0;
        end
        if ~isnan(partnerprob2(j))
            if isempty(find(cclust2(j,cclustID.clust)==cell2mat(stblclust2(1,:)),1))
                stblclust2{1,cclust2(j,cclustID.clust)} = cclust2(j,cclustID.clust);
                if partnerprob2(j)>0 
                    stblclust2{2,cclust2(j,cclustID.clust)} = cclust1(j,cclustID.clust);
                end
                stblclust2{3,cclust2(j,cclustID.clust)} = length(find(cclust2(j,cclustID.clust)==cclust2(:,cclustID.clust)));
                stblclust2{4,cclust2(j,cclustID.clust)} = 1;
            elseif partnerprob2(j)>0               
                stblclust2{2,cclust2(j,cclustID.clust)} = [stblclust2{2,cclust2(j,cclustID.clust)} cclust1(j,cclustID.clust)];
                stblclust2{4,cclust2(j,cclustID.clust)} = stblclust2{4,cclust2(j,cclustID.clust)} +1;
            end  
%         else
%             stblclust2{1,j} = j;
%             stblclust2{2,j} = 0;
%             stblclust2{3,j} = 0;
%             stblclust2{4,j} = 0;
        end        
    end
    partnerhist2(i,:) = histcounts(partnerprob2,0:.1:1,'normalization','probability'); %(partnerprob>0)
    
    for j = 1:size(stblclust2,2)
        stblclust2{5,j} = sum(isplace2(cclust2(:,cclustID.clust) == j));
        stblclust2{6,j} = sum(isstim2(cclust2(:,cclustID.clust) == j));
        stblclust2{7,j} = sum(isboth2(cclust2(:,cclustID.clust) == j));
    end
    
    %% mean r-value per day and cluster

    netcorr = CAIMcorr(expID(i+1),mouse).network.netcorr;
    netcorr(1:1+size(netcorr,1):end) = NaN;
    netcorr2 = CAIMcorr(trial,mouse).network.netcorr;
    netcorr2(1:1+size(netcorr2,1):end) = NaN;
    netcorr2(end+1:size(netcorr,1),end+1:size(netcorr,1)) = NaN;
%     e = 0;
%     A = [];
    for j = 1:size(stblclust2,2)
        %%
        k = max(cclust2(:,cclustID.clust))+1-j;
        b = cclust2(:,cclustID.clust) == k;
%         A(e+1:e+sum(b),e+1:e+sum(b)) = netcorr(b,b);
%         e = e+sum(b);
        stblclust2{8,k} = nanmean(nanmean(netcorr(b,b)));      
        stblclust2{9,k} = nanmean(nanmean(netcorr2(b,b)));
        c = zeros(1,size(stblclust2,2)-1);
        for kk = 1:size(stblclust2,2)
            if kk~=k
                bb = cclust2(:,cclustID.clust) == kk;
                c(kk) = nanmean(nanmean(netcorr2(b,bb))); 
            end
        end
        stblclust2{10,k} = c;        
    end
    
    %%
    clustsize = zeros(5,max([size(stblclust,2) size(stblclust2,2)]));
    clustcolor = [.5 .5 .5 ; 0 176/255 240/255; 1 0 1; 1 1 0];
    for j=1:max([size(stblclust,2) size(stblclust,2)])
        if j<=size(stblclust,2) && ~isempty(stblclust{1,j})
            clustsize(1,j) = cell2mat(stblclust(3,j));
            % color due to concept cells
            if stblclust{5,j} >0 && stblclust{6,j} == 0 && stblclust{7,j} == 0
                clustsize(5,j) = 1;
            elseif stblclust{5,j} ==0 && stblclust{6,j} > 0 && stblclust{7,j} == 0
                clustsize(5,j) = 2;
            elseif (stblclust{5,j} > 0 && stblclust{6,j} > 0) || stblclust{7,j} > 0
                clustsize(5,j) = 3;
            end
        else
            clustsize(1,j) = NaN;
        end
        if j<=size(stblclust2,2) && ~isempty(stblclust2{1,j})
            clustsize(2,j) = cell2mat(stblclust2(3,j));
            % color due to concept cells
            if stblclust2{5,j} >0 && stblclust2{6,j} == 0 && stblclust2{7,j} == 0
                clustsize(6,j) = 1;
            elseif stblclust2{5,j} ==0 && stblclust2{6,j} > 0 && stblclust2{7,j} == 0
                clustsize(6,j) = 2;
            elseif (stblclust2{5,j} > 0 && stblclust2{6,j} > 0) || stblclust2{7,j} > 0
                clustsize(6,j) = 3;
            end
        else
            clustsize(2,j) = NaN;
        end
        clustsize(3,j) = nansum(clustsize(1,1:j-1));
        clustsize(4,j) = nansum(clustsize(2,1:j-1));
        
    end

    figure('color',[1 1 1],...
        'position',[500 50 1.5*[420 594]],...
        'renderer','painters')
    hold on
    %bar(clustsize(1:2,:),.1,'w','stacked');
    
    % Moving partners
%     for j = 1:size(stblclust,2)     
%         b = stblclust{2,j};
%         if length(b)>1
%             b = sort(b);
%             for k = 1: length(b)                             
%                 aa = clustsize(3,j);
%                 bb = clustsize(4,b(k));
%                 plot([1+clustsize(1,j)/2/100 2-clustsize(2,b(k))/2/100],[aa bb],'linewidth',1,'color',[.6 .6 .6]) 
%                 clustsize(3,j) = aa + 1;
%                 clustsize(4,b(k)) = bb + 1;                
%             end
%         end
%     end
    
    % Moving concept cells
%     isthere = ~isnan(cclust1(:,cclustID.clust)) & ~isnan(cclust2(:,cclustID.clust));
%     isplace = [isthere & isplace1   isthere & isplace2];
%     isstim = [isthere & isstim1   isthere & isstim2];
%     isboth = [isthere & isboth1   isthere & isboth2];
%     
%     for j = 1:size(cclust1,1)
%         if isthere(j) 
%             a = cclust1(j,cclustID.clust);
%             aa = clustsize(3,a);
%             b = cclust2(j,cclustID.clust);
%             bb = clustsize(4,b);
%             if isplace(j,1)    
%                 plot([1+clustsize(1,a)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(2,:))                         
%             elseif isstim(j,1) 
%                 plot([1+clustsize(1,a)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(3,:)) 
%             elseif isboth(j,1) 
%                 plot([1+clustsize(1,a)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(4,:))  
%             elseif isplace(j,2) || isstim(j,2) || isboth(j,2)
%                 plot([1+clustsize(1,a)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(1,:))
%             end
%             if isplace(j,2)
%                 plot([1.5 2-clustsize(2,b)/2/100],[(bb+aa)/2 bb],'color',clustcolor(2,:))
%             elseif  isstim(j,2)                         
%                 plot([1.5 2-clustsize(2,b)/2/100],[(bb+aa)/2 bb],'color',clustcolor(3,:))
%             elseif  isboth(j,2)
%                 plot([1.5 2-clustsize(2,b)/2/100],[(bb+aa)/2 bb],'color',clustcolor(4,:))
%             elseif isplace(j,1) || isstim(j,1) || isboth(j,1)
%                 plot([1.5 2-clustsize(2,b)/2/100],[(bb+aa)/2 bb],'color',clustcolor(1,:))
%             end 
%             if isplace(j,2) || isstim(j,2) || isboth(j,2)
%                 clustsize(3,a) = aa + 1;
%                 clustsize(4,b) = bb + 1;
%             end
%         end
%     end

    isthere = ~isnan(cclust1(:,cclustID.clust)) & ~isnan(cclust2(:,cclustID.clust));
    isplace = [isthere & isplace1   isthere & isplace2];
    isstim = [isthere & isstim1   isthere & isstim2];
    isboth = [isthere & isboth1   isthere & isboth2];

    for j = 1:size(stblclust,2)
        a = find(cclust1(:,cclustID.clust)==j);
        b = zeros(length(a),7);
        for k =1:length(a)
            if isplace(a(k),1)
                b(k,1) = cclust2(a(k),cclustID.clust);
                b(k,2) = 1;
            elseif isstim(a(k),1)             
                b(k,1) = cclust2(a(k),cclustID.clust);
                b(k,4) = 1;
            elseif isboth(a(k),1) 
                b(k,1) = cclust2(a(k),cclustID.clust);
                b(k,6) = 1; 
            elseif ~isempty(find(cclust2(a(k),cclustID.clust)==stblclust{2,j},1))
                b(k,1) = cclust2(a(k),cclustID.clust);
            end
            if isplace(a(k),2)
                b(k,1) = cclust2(a(k),cclustID.clust);
                b(k,3) = 1;
            elseif isstim(a(k),2)             
                b(k,1) = cclust2(a(k),cclustID.clust);
                b(k,5) = 1;
            elseif isboth(a(k),2) 
                b(k,1) = cclust2(a(k),cclustID.clust);
                b(k,7) = 1; 
            elseif ~isempty(find(cclust2(a(k),cclustID.clust)==stblclust{2,j},1))
                b(k,1) = cclust2(a(k),cclustID.clust);
            end
        end
        b = b(sum(b,2)>0,:);
        [~,c] = sort(b(:,1));
        b = b(c,:);
        
    %%
        for k = 1 : size(b,1)                             
            aa = clustsize(3,j);
            bb = clustsize(4,b(k,1));
            if b(k,2)    
                plot([1+clustsize(1,j)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(2,:))                         
            elseif b(k,4) 
                plot([1+clustsize(1,j)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(3,:)) 
            elseif b(k,6) 
                plot([1+clustsize(1,j)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(4,:))  
            else
                plot([1+clustsize(1,j)/2/100 1.5],[aa (bb+aa)/2],'color',clustcolor(1,:))
            end
            if b(k,3)
                plot([1.5 2-clustsize(2,b(k,1))/2/100],[(bb+aa)/2 bb],'color',clustcolor(2,:))
            elseif  b(k,5)                         
                plot([1.5 2-clustsize(2,b(k,1))/2/100],[(bb+aa)/2 bb],'color',clustcolor(3,:))
            elseif  b(k,7)
                plot([1.5 2-clustsize(2,b(k,1))/2/100],[(bb+aa)/2 bb],'color',clustcolor(4,:))
            else
                plot([1.5 2-clustsize(2,b(k,1))/2/100],[(bb+aa)/2 bb],'color',clustcolor(1,:))
            end        
            clustsize(3,j) = aa + 1;
            clustsize(4,b(k,1)) = bb + 1;           
        end
    end
    % Cluster size and concepts
    for j = 1:size(stblclust,2)       
        if j<=size(stblclust,2)
            barwidth = clustsize(1,j)/100;
            fill([1-barwidth/2 1+barwidth/2 1+barwidth/2 1-barwidth/2],[sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j)) sum(clustsize(1,1:j))],clustcolor(1,:))
            b = 0;
            if stblclust{5,j}>0
                a = barwidth * stblclust{5,j}/stblclust{3,j};
                a = [1-barwidth/2+b 1-barwidth/2+a+b ];
                b = barwidth * stblclust{5,j}/stblclust{3,j};
                fill([a(1) a(2) a(2) a(1)],[sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j)) sum(clustsize(1,1:j))],clustcolor(2,:))
            end
            if stblclust{6,j}>0
                a = barwidth * stblclust{6,j}/stblclust{3,j};
                a = [1-barwidth/2+b 1-barwidth/2+a+b ];
                b = barwidth * stblclust{6,j} /stblclust{3,j};
                fill([a(1) a(2) a(2) a(1)],[sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j)) sum(clustsize(1,1:j))],clustcolor(3,:))
            end
            if stblclust{7,j}>0
                a = barwidth * stblclust{7,j}/stblclust{3,j};
                a = [1-barwidth/2+b 1-barwidth/2+a+b ];
                fill([a(1) a(2) a(2) a(1)],[sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j)) sum(clustsize(1,1:j))],clustcolor(4,:))
            end
%             plot([1 1],[sum(clustsize(1,1:j-1)) sum(clustsize(1,1:j))],'linewidth',10,'color',clustcolor(clustsize(5,j)+1,:))      
        
        end
        if j<=size(stblclust2,2)
            barwidth = clustsize(2,j)/100;
            fill([2-barwidth/2 2+barwidth/2 2+barwidth/2 2-barwidth/2],[sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j)) sum(clustsize(2,1:j))],clustcolor(1,:))
            b = 0;
            if stblclust2{5,j}>0
                a = barwidth * stblclust2{5,j}/stblclust2{3,j};
                a = [2-barwidth/2+b 2-barwidth/2+a+b ];
                b = stblclust2{5,j}/stblclust2{3,j}*barwidth;
                fill([a(1) a(2) a(2) a(1)],[sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j)) sum(clustsize(2,1:j))],clustcolor(2,:))
            end
            if stblclust2{6,j}>0
                a = barwidth * stblclust2{6,j}/stblclust2{3,j};
                a = [2-barwidth/2+b 2-barwidth/2+a+b ];
                b = stblclust2{6,j}/stblclust2{3,j}*barwidth;
                fill([a(1) a(2) a(2) a(1)],[sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j)) sum(clustsize(2,1:j))],clustcolor(3,:))
            end
            if stblclust2{7,j}>0
                a = barwidth * stblclust2{7,j}/stblclust2{3,j};
                a = [2-barwidth/2+b 2-barwidth/2+a+b ];
                fill([a(1) a(2) a(2) a(1)],[sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j)) sum(clustsize(2,1:j))],clustcolor(4,:))
            end
%             plot([2.2 2.2],[sum(clustsize(2,1:j-1)) sum(clustsize(2,1:j))],'linewidth',10,'color',clustcolor(clustsize(6,j)+1,:))      
        end
    end
    axis off
    hold off
    mouseID = CAIMcorr(trial,mouse).cclust(:,1);
    mouseID = mouseID(~isnan(mouseID));
    mouseID = num2str(mouseID(1));
    title(['M'  mouseID ', ' experiment{trial} ' to ' experiment{expID(i+1)}])
%     printpdf(['M' mouseID ',' experiment{trial}])
  %%  
    ClustProp = [ClustProp stblclust2];     
end

%% Cluster Correlation plots with sorting over days

% cclust = cell(length(expID),1);
% netraster = cell(length(expID),1);
% netnum = zeros(length(expID)+1,1);
% isplace = cell(length(expID),1);
% isstim = cell(length(expID),1);
% isboth = cell(length(expID),1);
% isthere = cell(length(expID),1);
% 
% for i = 1:length(expID)
%     trial = expID(i);
%     netraster{i} = CAIMcorr(trial,mouse).network.netraster;
%     netnum(i+1) = size(netraster{i},2)+netnum(i);
%     cclusttemp = CAIMcorr(trial,mouse).cclust;
%     
% %     active field sorting
%     netraster{i}(CAIMcorr(trial,mouse).inField==0,:) = 0;
%     cclusttemp(CAIMcorr(trial,mouse).inField==0,:) = NaN;
% 
%     isplace{i} = (cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & ~(cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
%     isstim{i}  = ~(cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & (cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
%     isboth{i}  = (cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & (cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
%     isthere{i} = ~isnan(cclusttemp(:,1));
%     cclust{i} = cclusttemp;
% end
% 
% for j=1:length(expID)-1
%     %%
%     netraster1 = netraster{j};
%     netraster2 = netraster{j+1};
%     netraster2 = netraster2(1:size(netraster1,1),:);
%     a = sum(netraster2,2)>0;
%     netraster1 = netraster1(a,:);
%     netraster2 = netraster2(a,:);
%     figure
%     
%     subplot(1,2,1)
%     [a,c] = netcorrplot(netraster1,[],[]);
%     title(experiment(expID(j)))
%     subplot(1,2,2)
%     netcorrplot(netraster2,a,c);
%     title(experiment(expID(j+1)))
% end


function [a,c] = netcorrplot(netraster,a,c)
    
    if isempty(c)
        [netcorr,~] = corr(netraster');
        a = ~isnan(diag(netcorr));
        netcorr = netcorr(a,a);
        Y = netcorr;
        Z = linkage(Y,'weighted','seuclidean');
        
%         D = pdist(Y,'seuclidean');     
%         c = optimalleaforder(Z,D);
        
        [aa,bb] = histcounts(nansum(netraster),'normalization','probability');
        numcell = sum(a);
        sumaa = 0;
        for ii = 1:length(aa)
            sumaa = sumaa + aa(ii);
            if sumaa >0.5
                meannet = bb(ii)+.5;
                break
            end
        end 
        numclust = round(numcell/meannet);
        C = cluster(Z,'maxclust',numclust);
        [~,c] = sort(C,'descend');
               
    else
        [netcorr,~] = corr(netraster');
        netcorr = netcorr(a,a);      
        Y = netcorr;      
    end
    
    A = Y(c,c);
    imagesc(A)
    colormap(jet)

end
end