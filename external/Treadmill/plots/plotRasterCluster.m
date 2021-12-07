% raster scatter color coded events
experiment = {'Baseline'};
mouse = {'177', '184', '446', '453'}; 
%load('/media/2Photon/Nicola/AnalysisFINAL/BigFatSummary.mat')
load('cclustID.mat')
%%
plotk = 11;
num = [3:5 8];
% x =scn.tsscn/1000/60;
for iii = num
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(iii,fullCAIM(j)).A)
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
    fullCAIM(emptyCAIM) = [];
   
    for jjj = 1:length(fullCAIM)    
        kkk = fullCAIM(jjj);
        
        

        % S = logical(caim.S_bin);
        S = CAIM(iii,kkk).S;
        % PcID = scn.cellID;
        % PcID = PcID(PcID(:,6)<.05,1);
        PcID = CAIM(iii,kkk).cclust;
        PcID = find(PcID(:,cclustID.plcvct)>0 & PcID(:,cclustID.plcvctp)<.05);
     
        tsscn = CAIM(iii,kkk).behave.tsscn;
        running = CAIM(iii,kkk).behave.running;
        network = CAIM(iii,kkk).network;
        speed = CAIM(iii,kkk).behave.speed; 
        % network = caim.network;
        time = [1 network.netpos(end)];
        j = 1;
        netraster = zeros(size(network.netraster,1),1);
        netpos = zeros(1,1);
        if network.cellID ~= 0
            c = network.cellID(:,2);
            c(isnan(c)) = max(c)+1:length(c);
            clust = network.cellID(:,1);
            clust(isnan(clust)) = max(clust)+1;
            [~,c] = sort(clust,'ascend');
            % c = 1:length(c);
            for i = 1:length(network.netID)
                if network.netpos(i)>time(1) && network.netpos(i) < time(2) 
                    netpos(j) = network.netpos(i);
                    netraster(:,j) = network.netraster(c,i);
                    j = j+1;
                end
            end

            % for i = 1:2:size(netraster,1)  
            %     plot([0 size(netraster,2)+1],[i i] ,':','color',[.7 .7 .7])
            % end
            figure('color',[1 1 1],...
                'renderer','painters',...
                'visible','on',...
                'Units','centimeters',...
                'position',[10 2 [ 2.5*sqrt(2)*8.9 2.5*8.9]],...
                'PaperUnits','centimeters',...
                'PaperSize', [2.5*sqrt(2)*8.9 2.5*8.9])

            netcolor = prism(max(clust)+1);
            cellim = 160;if cellim>size(netraster,2);cellim = size(netraster,2);end
            j = 1;
            for i = 1:cellim% 
                plot([netpos(i)/1000 netpos(i)/1000],[0 500] ,':','color',[.7 .7 .7])
                hold on
            end

            clustamp = zeros(max(clust),cellim);
            for ii = 1:max(clust)
                for i = 1:cellim% 
                    k = find(netraster(:,i) >= 1 & clust(c)==ii);        
                    clustamp(ii,i) = length(k);
                    scatter(ones(1,length(k))*netpos(i)/1000,k,10,'filled','MarkerFaceColor',netcolor(ii,:));       
                end
            end

            ax = gca;
            ax.FontSize = ftsz-2;
            ax.LineWidth = lnwd;
            ax.YColor = [0 0 0];
            ax.XColor = [0 0 0];
            ax.Color = [1 1 1];
            ax.YTick = [0:100:500];
            ax.XLim = [time(1) time(2)]/1000;
            ax.YLim = [1 size(netraster,1)];
            ylabel('cell ID','Fontsize',ftsz-2)
            xlabel('time/s','Fontsize',ftsz-2)

            for i = 1:length(PcID)
            %     pS = scn.tsscn(S(PcID(i),:)==1 & scn.running'==1)/1000;
                pS = tsscn(S(PcID(i),:)==1 & running'==1)/1000;
                scatter(pS,c(PcID(i))*ones(length(pS),1),10,'filled','MarkerFaceColor',[0 0 0])
            end
            hold off
            title([mouse{kkk} ' ' experiment{iii}])
            %%
            plotk = plotk+1;
            printpdf([num2str(plotk)])
        end
    end
end

files = dir('*.pdf');
append_pdfs(['Raster plots.pdf'],files.name)
delete(files.name)
close all
%%
smoothamp = zeros(size(clustamp));
for i = 1:size(clustamp,1)
    smoothamp(i,:) = clustamp(i,:);
    smoothamp(i,:) = smooth(smoothamp(i,:),5);
    smoothamp(i,:) = mat2gray(smoothamp(i,:));
    
end
figure
imagesc(smoothamp)
colormap(jet)