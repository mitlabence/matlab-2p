% Examplefigures

% pathname1 = 'Z:\imaging Negar & Martin\M103\Analysis\Baseline\';
% pathname1 = 'Z:\imaging Negar & Martin\M103\Analysis\Airpuff\';
pathname1 = 'Z:\imaging Negar & Martin\M226\Analysis\Airpuff\';
[files1,samedate] = FindDataSets(pathname1);
[belt1,caim1] = ConnecCaim(pathname1,files1,samedate{2}); % 5 for baseline
[belt,scn1] = placecor(belt1,caim1,[],10);
caim1.network = netevents(caim1,scn1);
caim1.fireprop = firefreq(caim1,scn1);
scn1 = stimcor(belt,caim1,scn1);
caim1.network = ClustAna(caim1,scn1);

%%
pathname = 'Z:\imaging Negar & Martin\M234a\Analysis\Baseline\';
trial = 2;
[~,~,prefiles] = FindDataSets(pathname);
disp(['load ' prefiles(trial).name])
load(prefiles(trial).name)
%%
lnwd = 3;
ftsz = 25;

%% Network overview picture
figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   350   420])
hold on

netcolor = [1 0 0;1 0 1; 1 1 0;0 1 0;0 0 .8;.5 0 1;1 .5 0;0 1 1;1 0 .5;0 .5 0;.8 .5 0;0 .8 .5;];

time = [4800 6000];
b = zeros(caim1.options.d1,caim1.options.d2);
j = 1;
for i = 1:length(caim1.network.netID)
    if caim1.network.netpos(i)>time(1) && caim1.network.netpos(i) < time(2)    
        a = full(reshape(sum(caim1.A(:,caim1.network.netID{i}),2),caim1.options.d1,caim1.options.d2));
        a = a(:,:,[1 1 1])/max(max(a));       
        a(:,:,netcolor(j,:)==0) = 0;
        j = j+1;      
        if  j == length(netcolor)
            j = 1;
        end  
        b = b + a;  
    end  
end
%     imshow(1.1*Cn+1*b);

C = imfuse(1.2*caim1.Cn+.6,.8*b,'blend','Scaling','joint');
% a = sum(b,3)==0;a = double(a);
% b(:,:,1) = b(:,:,1)+a;
% b(:,:,2) = b(:,:,2)+a;
% b(:,:,3) = b(:,:,3)+a;

% C = C(100:300,100:450,:);
imshow(1.8*C)

% imagesc(C)
% ax = gca;
% ax.XLim = [100 400];%[0 size(netraster,2)];
% ax.YLim = [100 250];
% axis off

hold off
exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\Network ID')
%% rasterplot example
figure('name','Rasterplo',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   650   650])
hold on
j = 1;
netraster = zeros(size(caim1.network.netraster,1),1);
netpos = zeros(1,1);
for i = 1:length(caim1.network.netID)
    if caim1.network.netpos(i)>time(1) && caim1.network.netpos(i) < time(2) 
        netpos(j) = caim1.network.netpos(i);
        netraster(:,j) = caim1.network.netraster(:,i);
        j = j+1;
    end
end

% for i = 1:2:size(netraster,1)  
%     plot([0 size(netraster,2)+1],[i i] ,':','color',[.7 .7 .7])
% end

cellim = 88;if cellim>size(netraster,2);cellim = size(netraster,2);end
j = 1;
for i = 1:cellim%   
    k = (find(sum(netraster(:,i,:),3) >= 1));
    scatter(ones(1,length(k))*scn1.tsscn(netpos(i))/1000,k,30,'filled','MarkerFaceColor',netcolor(j,:));
    j = j +1;
    if j > size(netcolor,1)
        j = 1;
    end
end

ax = gca;
ax.FontSize = ftsz-10;
ax.XLim = [320 400];%[0 size(netraster,2)];%
ax.YLim = [0 size(netraster,1)];
ax.LineWidth = lnwd;
% ax.YTick = 0:10:size(netraster,1);
ax.XTick = ax.XLim(1):10:ax.XLim(2);
ax.XTickLabel = 0:10:ax.XLim(2)-ax.XLim(1);
% ax.YTickLabel = {};
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
xlabel('time /s','FontSize',ftsz-2)
ylabel('cell #','FontSize',ftsz-2)
hold off
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\rasterplot')

%% Correlation matrix
figure('name','Netcorr',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   650   650])
a = caim1.network.netcorr;
a = a(~isnan(diag(a)),~isnan(diag(a)));
a(a<0) = 0;
imagesc(a)
colormap(jet)
ax = gca;
ax.FontSize = ftsz-10;

ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
xlabel('cell #','FontSize',ftsz-2)
ylabel('cell #','FontSize',ftsz-2)
hold off
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\rasterplot')
%% Cluster Matrix
figure('name','Netcorr',...
    'color',[1 1 1],...
    'renderer','painters',...
    'position',[550   150   650   650])
a = caim1.network.netcluster;
a = rot90(rot90(a));
% a = a(~isnan(diag(a)),~isnan(diag(a)));
a(a<0) = 0;
imagesc(a)
colormap(jet)
ax = gca;
ax.FontSize = ftsz-10;

ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];
xlabel('cell #','FontSize',ftsz-2)
ylabel('cell #','FontSize',ftsz-2)
hold off
% axis off
c = colorbar;
c.Color = [0 0 0];
c.FontSize = ftsz-10;
exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\cluster')

%% dendrogramm (to insert in ClustAna)
% ClustAna(caim1,scn1);
% figure('name','Netcorr',...
%     'color',[1 1 1],...
%     'renderer','painters',...
%     'position',[550   150   650   650])
% D = pdist(Y);
% cc = optimalleaforder(Z,D);
% dendrogram(Z,0,'Orientation','left','reorder',c,'ColorThreshold',0.9*max(Z(:,3)))
% ax = gca;
% ax.Color = [0 0 0];
% axis off
% % exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\dendrogram')
%% Cluster Matrix with features
figure('name','Netcorr',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   650   650])

a = caim1.network.netcluster;
% a = a(~isnan(diag(a)),~isnan(diag(a)));
a(a<0) = 0;
imagesc(a)
colormap(gray)
ax = gca;
ax.FontSize = ftsz-10;
ax.XLim = [200 500];
ax.YLim = ax.XLim;
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
xlabel('cell #','FontSize',ftsz-2)
ylabel('cell #','FontSize',ftsz-2)
hold off
% c = colorbar;
% c.Color = [1 1 1];
% c.FontSize = ftsz-10;

% if isfield(scn1,'plcfield') && ~isempty(scn1.plcfield)
%     b = caim1.network.cellID(:,2);
%     isplace = false(1,size(caim1.Y,1));
%     isplace(scn1.cellID(scn1.cellID(:,6)<.05,1)) = 1;
%     isplace = isplace(~isnan(b));
%     isplace = isplace(b(~isnan(b)));
%     isplace = find(isplace);
%     hold on
%     for i = isplace
%         scatter(i,i,25,'filled',...
%             'LineWidth',1,...
%             'MarkerEdgeColor',[0 1 1],...        
%             'MarkerFaceColor',[0 1 1])
%     end
%     if isfield(scn1,'airpuff')  
%         stim = scn1.airpuff;
%         isstim = false(1,size(caim1.Y,1));
%         isstim(stim.cellID(:,1)) = 1;
%         isstim = isstim(~isnan(b));
%         isstim = isstim(b(~isnan(b)));
%         isstim = find(isstim);
%         for i = isstim
%             scatter(i,i,25,'filled',...
%             'LineWidth',1,...
%             'MarkerEdgeColor',[1 0 1],...        
%             'MarkerFaceColor',[1 0 1])
%         end
%     end
%     hold off
% end
% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\cluster feature 2')
%% Polar plots of place cells

figure('name','Netcorr',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   650   850])
cells = [35 16 41];

% hold on
for j = 1:length(cells)
    subplot(3,2,2*j-1)
    i = find(scn.cellID(:,1)==cells(j));
    d = (max(scn.rounds)+1)*scn.spcpol{i,1};
    c = scn.spcpol{i,3};       
%     title(['Cell ' num2str(i)],...
%         'fontsize',7)
    polarplot(scn.theta,scn.rho,'color',[0 176/255,240/255]) 
    hold on               
    polarplot(c(:,1),c(:,2),'o','color',[0 1 0])
%     polarplot(d,'*');
    polarplot([0 real(-1j*log(d))],[0 abs(d)],...
        'LineWidth',lnwd,...
        'color',[1 0 1]);
    axis off
    hold off

    subplot(3,2,2*j)
    
    a = length(find(scn.shufpol(i,2,:)>scn.spcpol{i,2}))/(size(scn.shufpol(i,2,:),3));
    histogram(scn.shufpol(i,2,:),...
        'Normalization','Probability',...
        'FaceColor',[.9 .9 .9])
    box off
    hold on    
    plot([scn.spcpol{i,2} scn.spcpol{i,2}],[0 .2],...
        'LineWidth',lnwd,...
        'color',[0 1 1])
    xlim([0 1])
    ax = gca;
    ax.FontSize = ftsz-5;
    ax.YColor = [1 1 1];
    ax.XColor = [1 1 1];
    ax.Color = [0 0 0];
    ax.YLim = [0 .2];
    ax.LineWidth = lnwd;
    hold off
    if scn.cellID(i,6) == 0
        text(.7,.17,'p < e^-^6','FontSize',ftsz-5,'color',[1 1 1])
    else
        text(.7,.17,['p = ' num2str(round(100*scn.cellID(i,6))/100) ],'FontSize',ftsz-5,'color',[1 1 1])
    end

end

% exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\place')
%%
for i = 1:length(cells)
    figure('name','Netcorr',...
        'color',[0 0 0],...
        'renderer','painters',...
        'position',[550   150   650   80])

    j = find(scn.cellID(:,1)==cells(i));
    imagesc(scn.plcfield(j,:,1))
    axis off
%     exportfig(['C:\Users\Admin\Dropbox\Lab\Vorträge\place ' num2str(i)])
end   

%%
pathname = 'Z:\imaging Negar & Martin\M234\Analysis\Airpuff\';
trial = 2;
[~,~,prefiles] = FindDataSets(pathname);
disp(['load ' prefiles(trial).name])
load(prefiles(trial).name)

%% Example AP resposnsive cell M234.1 cell 2 M234.2 cell 9 

stim = scn.airpuff;
i = 9;% CellID
% subplot(10,3,[5 8 11 14])
axes('Position',[.41 .54 .215 .3])
for j = 1:size(stim.resp,2)
    scatter(stim.times(j,stim.resp(i,j,:)==1)/1000,stim.resp(i,j,stim.resp(i,j,:)==1)+j-1,'g','filled') 
    hold on
end


ax = gca;
ax.XLim = [-1 3];
ax.YLim = [0 size(stim.resp,2)];
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.Color = [1 1 1];     
ax.LineWidth = lnwd;
ax.FontSize = ftsz-5;
box('off')
plot([0 0],[0 size(stim.resp,2)],'m','linewidth',lnwd)
hold off
ylabel('Airpuff #','FontSize',ftsz)

title('b)','Units','centimeter','position',[-1 6.8],'FontSize',ftsz)

% subplot(10,3,2)
axes('Position',[.41 .875 .215 .1])
imagesc(stim.times(i,:)./1000,1,permute(sum(stim.resp(i,:,:)),[1 3 2]))
hold on
ax = gca;
ax.XLim = [-1 3];
ax.YLim = [0 1];
% ax.YColor = [0 0 0];
% ax.XColor = [0 0 0];
% ax.Color = [1 1 1];     
% ax.LineWidth = lnwd;
% ax.FontSize = ftsz-5;
box('off')
% plot([0 0],[0 1],'m','linewidth',lnwd)
hold off
axis off
% ylabel('Airpuff #','FontSize',ftsz)

%%
figure('name','Firing Frequency',...
    'color',[0 0 0],...
    'renderer','painters',...
    'position',[550   150   520   620])

stim = scn.airpuff;
% stim = scn.network;
cclusttemp = CAIM(9,1).cclust;
a = zeros(size(stim.resp,1),size(stim.resp,3)); 
pvalues = zeros(1,size(stim.resp,1));

for j = 1:size(stim.resp,1)        
    pvalues(j) =  cclusttemp(stim.cellID(j,1),cclustID.airpp);
    a(j,:) = sum(stim.resp(j,:,:),2);
%     a(j,:) = (a(j,:)-min(a(j,:)))./max(a(j,:)-min(a(j,:)));
%     a(j,:) = smooth(a(j,:));
end

aa = a(pvalues>.05,:);
a = a(pvalues<=.05,:);


subplot(8,1,1:2)
imagesc(stim.times(1,:)/1000,1:size(a,1),a)
box('off')
%%
hold on 
plot([0 0],[0 size(a,1)+.5],'m','linewidth',lnwd)

ax = gca;
ax.XLim = [-1 3];
% ax.YLim = [0 size(stim.resp,1)+20];
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];
ax.XTick = [];
ax.LineWidth = lnwd-1;
ax.FontSize = ftsz-10;
% ylabel('Airpuff #','FontSize',ftsz)

hold off

stim = scn.airpuff;

a = zeros(size(stim.nr.resp,1),size(stim.nr.resp,3)); 
for j = 1:size(stim.nr.resp,1)        
    a(j,:) = sum(stim.nr.resp(j,:,:),2);
%     a(j,:) = (a(j,:)-min(a(j,:)))./max(a(j,:)-min(a(j,:)));
%     a(j,:) = smooth(a(j,:));
end        

a = [aa ;a];
subplot(8,1,3:8)

imagesc(stim.times(1,:)/1000,1:size(a,1),a)
hold on 
plot([0 0],[0 size(a,1)+.5],'m','linewidth',lnwd)
box('off')

xlabel('time /s','FontSize',ftsz)
ylabel('Cell ID','FontSize',ftsz)
ax = gca;
ax.XLim = [-1 3];
% ax.YLim = [0 200];%size(stim.resp,1)+20];
ax.YColor = [1 1 1];
ax.XColor = [1 1 1];
ax.Color = [0 0 0];     
ax.LineWidth = lnwd-1;
ax.FontSize = ftsz-10;


hold off
%  exportfig('C:\Users\Admin\Dropbox\Lab\Vorträge\AP trial example')