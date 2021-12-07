ftsz = 8;
lnwd = 1;
numice = [1:11];
% The rows refer to experimental conditions the were pooled. Decide here
% which experiments you want to pool
numexp = [5];

plcfieldD = [];
plccenterD = [];
plcfieldV = [];
plccenterV = [];
            
% This loop goes through the rows (the experiemts)
for i = 1:length(numexp)
    
    % first sort out mice that might not contain the desired variable or
    % are excluded
    fullCAIM = numice;
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if ~isfield(CAIM(numexp(i),fullCAIM(j)),'cclust')
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    fullCAIM(emptyCAIM) = [];    
    
    % Here happens the read out
    % The loop goes through the included mice
    
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j);  
        
        cclust = CAIM(numexp(i),k).cclust;
        isplace = cclust(:,cclustID.plcvct)>0 & cclust(:,cclustID.plcvctp)<=.05;
        
        plcfieldtemp = CAIM(numexp(i),k).plcfield(:,1:150);
        plcfieldtemp = plcfieldtemp(isplace,:);
        plccentertemp = cclust(isplace,cclustID.plcfld);
        
        if j < 6
            plcfieldD = [plcfieldD; plcfieldtemp];
            plccenterD = [plccenterD; plccentertemp];
        else
            plcfieldV = [plcfieldV; plcfieldtemp];
            plccenterV = [plccenterV; plccentertemp];
        end
    end
    
end
%%

figure
% sort angles

[plccenterD,aa] = sort(plccenterD);
plcfieldD = plcfieldD(aa,:);
c = find(abs(plccenterD)==min(abs(plccenterD)));
plcfieldD = plcfieldD([c:end 1:c-1],:);

imagesc(plcfieldD)
% pcolor(plcfield1); colormap(jet); shading interp
% cyan colormap

colormap(jet)

% colorbar
xlabel('position on belt /cm','Fontsize',ftsz-2)
ylabel('cell ID','Fontsize',ftsz-2)
ax = gca;
ax.XTick = [100/3 100*2/3 100];
ax.XTickLabel = [50 100 150];
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
% ax.XLim = [0 30];
% ax.YLim = [0 1];

%%
figure
% sort angles

[plccenterV,aa] = sort(plccenterV);
plcfieldV = plcfieldV(aa,:);
c = find(abs(plccenterV)==min(abs(plccenterV)));
plcfieldV = plcfieldV([c:end 1:c-1],:);

imagesc(plcfieldV)
% pcolor(plcfield1); colormap(jet); shading interp
% cyan colormap

colormap(jet)

% colorbar
xlabel('position on belt /cm','Fontsize',ftsz-2)
ylabel('cell ID','Fontsize',ftsz-2)
ax = gca;
ax.XTick = [100/3 50*2/3 100];
ax.XTickLabel = [50 100 150];
ax.FontSize = ftsz-2;
ax.LineWidth = lnwd;
% ax.XLim = [0 30];
% ax.YLim = [0 1];