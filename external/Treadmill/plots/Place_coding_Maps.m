%%Place coding maps
numice = [1:11];
numexp = [5];

ftsz = 8;
figSize =[5 5 ];

figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[20 5 figSize ],...
    'PaperUnits','centimeters',...
    'PaperSize', figSize )

%use = [1]; % exp trials to compare
% Mice to use cells from
%mice = [103 155 158 194 195 224 226 227];


cellID = cclust(:,cclustID.expID,:);
cellID = max(cellID,[],3);
cellID = cellID == mice;
cellID = max(cellID,[],2);
thresh = 1;% threshold of appearance

samecelltemp = samecell(use,cellID);
actcell = samecelltemp;
actcell(actcell>0) = 1;
actcell = sum(actcell,1);

cclusttemp = cclust;
exclude = samecell' == 0;
for i = 1:size(cclusttemp,3)
    cclusttemp(exclude(:,i),:,i) = NaN;
end
cclusttemp = cclust(cellID,:,:);
cclusttemp = cclusttemp(actcell == 1,:,use);
samecelltemp = samecelltemp(:,actcell == 1);


a = cclusttemp(:,cclustID.plcvct,:)>0 & cclusttemp(:,cclustID.plcvctp,:)<=.05;
b = cclusttemp(max(a(:,1,:),[],3),:,:);

plcfieldtemp = plcfield(cellID,:,:);
plcfieldtemp = plcfieldtemp(actcell == 1,:,use);
plcfieldtemp = plcfieldtemp(max(a(:,1,:),[],3),:,:); 
for i = 1:size(plcfieldtemp,1)
    for j = 1:size(plcfieldtemp,3)
     plcfieldtemp(i,:,j) = (plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)))./max(plcfieldtemp(i,:,j)-min(plcfieldtemp(i,:,j)));
    end
end

% plcfields in basline condition
% sort angles
plccenter = b(:,cclustID.plcvctang,1);%plcfld
plccenter = plccenter(b(:,cclustID.plcvct,1)>0);
[plccenter,aa] = sort(plccenter);
plcfield1 = plcfieldtemp(b(:,cclustID.plcvct,1)>0,:,1);
plcfield1 = plcfield1(aa,:);
c = find(abs(plccenter)==min(abs(plccenter)))-7;
plcfield1 = plcfield1([c:end 1:c-1],:);

imagesc(plcfield1)
% pcolor(plcfield1); colormap(jet); shading interp
% cyan colormap
clmplg = floor(256*[0 1]);
mycolormap = zeros(256,3);
j = clmplg(2);

for i = 1:sum(clmplg)
    j = j-1;
    mycolormap(i,:) = [1-(clmplg(2)-j)/clmplg(2) 1 1];
end

% colormap(plotb1,mycolormap)
colormap(plotd1,jet)

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

