function PlotComponents(caim,nums,save,mov)

A = caim.A;
C = caim.C;
Cn = caim.Cn;
Df = caim.Df;
options = caim.options;
%%
tsscn = (1:length(C(1,:)))/15;
dt = tsscn(end)/length(tsscn);

if isempty(mov)
    scrsz = get(0,'ScreenSize');
    h = figure('Position',[scrsz(1) scrsz(4)/4 scrsz(3)-scrsz(1) 2*scrsz(4)/3],...
        'color',[1 1 1]);
    mov = length(C(1,:));
    pltcol = [0 0 0];
else
    h = gcf;
    set(h,'color',[0 0 0]);
    pltcol = [1 1 1];
end


ftsz = 20;

if isempty(nums)
     nums = 1:length(Df);
end

%% Plot traces
subplot(1,2,2)

for i = 1:length(nums)

%     aa = Y(nums(i),:)/Df(nums(i));
% %     aa = aa-min(aa);
%     plot(tsscn,aa+((i-1)),...
%         'color',[0 0 1],...
%         'linewidth',2);
    
    aa = C(nums(i),:)/Df(nums(i));
    aa = aa-min(aa);
    plot(tsscn(1:end),aa(1:end)+((i-1)),...
        'color',pltcol,...
        'linewidth',2);
    hold on
%     aa = S(nums(i),:)/max(S(nums(i),:));
%     plot(tsscn,aa+((i-1)),...
%         'color',[1 0 0],...
%         'linewidth',2);
    if mov < length(C(1,:))
        plot([tsscn(mov) tsscn(mov)],[0 length(nums)],...
            'color',[0 1 0],...
            'linewidth',2);
    end
end

set(gca,'ylim', [0 length(nums)]);
set(gca,'xlim', [0 tsscn(end)]);
axis off

%% plot stimulus (hard coded..sorry)

    plot([tsscn(1530) tsscn(1530)],[0 length(nums)],...
            'color',[0 0 1],...
            'linewidth',2);


%% scale bars
    fig_pos=get(gca,'position');
    yp1=fig_pos(2)-0.01;
    yp2=fig_pos(2)+fig_pos(4)/length(nums)-.01;%/(abs(min(aa))+max(aa));%+1/(length(nums)*max(a));
    xp1=fig_pos(1)-0.02;
    xp2=fig_pos(1)+fig_pos(3)/length(aa)/dt*10-0.02;


    h1 = annotation('textarrow',[xp1,xp2],[yp1,yp1],...
           'color',pltcol,...
           'HeadStyle','none',...       'HeadWidth',6,...       'Headlength',12,...
           'linewidth',2);

    h2 = annotation('textarrow',[xp1,xp1],[yp1+.18,yp1+.18],...
           'color',pltcol,...
           'HeadStyle','none',...
           'string','100% \DeltaF/F',...
           'fontsize',ftsz,...
           'VerticalAlignment','bottom',...  
           'TextRotation',90,...
           'HorizontalAlignment','right');

    h3 = annotation('textarrow',[xp1,xp1],[yp1,yp2],...
           'color',pltcol,...
           'HeadStyle','none',...
           'HeadWidth',6,...
           'Headlength',12,...
           'linewidth',2,...
           'string','             10 sec',...
           'fontsize',ftsz,...
           'VerticalAlignment','top',...
           'TextRotation',0,...
           'HorizontalAlignment','left');

    hold off
%% spatial components

    subplot(1,2,1)

    d1 = options.d1;
    d2 = options.d2;
    if max(Cn(:))>1
        imshow(Cn/max(max(Cn))); 
    else
        imshow(Cn)
    end
    hold on;
    
    if mov < length(C(1,:))
        a = zeros(d1,d2);
        for i = 1:length(A(1,:))  
            aa = full(reshape(A(:,i),d1,d2));
            a = a + aa(:,:,[1 1 1])/max(max(aa))*(C(i,mov)-min(C(i,:)))/(max(C(i,:)-min(C(i,:))));
        end
    else
        a = full(reshape(sum(A(:,nums),2),d1,d2));
%         a = full(reshape(sum(A(:,:),2),d1,d2));
        a = a(:,:,[1 1 1])/max(max(a));
    end
    
    a(:,:,[1 3]) = 0;
    faa = imshow(a);
    set(faa, 'AlphaData',a(:,:,2))
    
    pasplot = 0;
    if pasplot
        a = full(reshape(sum(caim.A_pas(:,:),2),d1,d2));
        a = a(:,:,[1 1 1])/max(max(a));
         a(:,:,[2 3]) = 0;
        faa = imshow(a);
        set(faa, 'AlphaData',a(:,:,1))
    end
    hold off

    %% scatter plot
    figure('position',[980 265 867 720],'color',[1 1 1])
    S = ones(size(caim.S_bin(nums,:)))-caim.S_bin(nums,:);
    imagesc(S)
    hold on
    colormap(gray)
    plot([1530 1530],[0 size(S,1)],...
        'color',[0 0 1],...
        'linewidth',2);
    
    set(gca,'XTickLabel',round((500:500:2500)*dt),...
        'fontsize',ftsz)
    hold off
    %% firing probabilty in bin

%     S = caim.S_bin(nums,:);

    bin = round(200*dt);
    numbin = floor(size(S,2)/bin);
    prob = zeros(size(S,1),numbin);
    for i = 1:numbin
        for j = 1:size(S,1)
            prob(j,i) = sum(S(j,(i-1)*bin+1:i*bin));
        end
    end
    % imagesc(prob)
    temp = (1:numbin)*bin*dt;


    figure('position',[680   529   947   449],'color',[1 1 1])

    % plot(temp,mean(prob),'color',[0 0 0])
    % scatter(temp,mean(prob),'b')
    errorbar(temp,mean(prob),std(prob)/length(nums),...
        'color',[0 0 0],...
        'linewidth',1.5)
    hold on
    plot([1530*dt 1530*dt],[0 ceil(max(mean(prob)))],...
            'color',[0 0 1],...
            'linewidth',2);
    grid on
    set(gca,'fontsize',ftsz,...
        'ylim', [0 ceil(max(mean(prob)))],...
        'xlim', [temp(1) temp(end)])
    hold off

    %% time to fire
    
    % nums = [1:9 12 14 15 19 20 21 24 25 26 28 34 51 65:71 74:81 83 84 88 89 90];
    
    firind = 1535;
    onset = zeros(1,size(S,1));
    for i = 1:size(S,1)
        if ~isempty(find(S(i,firind:end),1))
            onset(i) = find(S(i,firind:end),1)*dt;
        end
    end


    figure('color',[1 1 1],'position',[680   478   611   500])
    [counts,center] = hist(onset,20);
    cumset = cumsum(counts)/sum(counts);
    plot(center,cumset,'Linewidth',2)
    axis tight
    set(gca,'fontsize',ftsz,...
        'Linewidth',2)
    grid on

    figure('color',[1 1 1],'position',[680   465   169   513])
    boxplot(onset,'Colors','k',...
        'MedianStyle','line',...
        'width',1.5)

    set(gca,'fontsize',ftsz,...
        'Linewidth',2,...
        'XTick',[],'XTickLabel',[])
%%   
if save
    [filename, pathname] = uiputfile('*.png','Choose figure name','C:\Users\Admin\Desktop\');

    if ~filename == 0
        exportfig([pathname filename(1:end-4)])
    end
end
   

end
   