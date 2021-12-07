function PlotComponentsBehave(caim,nums,save,mov,Belt,ref,timeint,horizontal)

if nargin > 7 && length(horizontal) == 1 && horizontal == 0
    nxplot = 1;
    nyplot = 1;
    posplot = 0;
    ftsz = 12;
    lnwd = 1;
elseif length(horizontal) > 1
    nxplot = horizontal(1);
    nyplot = horizontal(2);
    posplot = horizontal(3);
    ftsz = 12;
    lnwd = 1;
else
    nxplot = 2;
    nyplot = 1;
    horizontal = true;  
    posplot = 0;
    ftsz = 12;
    lnwd = 1;
end

A = caim.A;
C = caim.C;
if isempty(ref)
    Cn = caim.Cn;
else
    Cn = ref;
end
Df = caim.Df;
% S = caim.S;
% Y = caim.Y;
options = caim.options;

if isempty(Belt)
    tsscn = (1:length(C(1,:)))/30;
else
    tsscn = Belt.tsscn;
end
dt = tsscn(end)/length(tsscn);

if isempty(timeint)
   timeint = [0 tsscn(end)];  
else
%    timeint = timeint*60*1000; % time interval to be plotted
end

if isempty(mov)
    if horizontal == 1      
        scrsz = get(0,'ScreenSize');
        h = figure('Position',[scrsz(1) scrsz(4)/4 scrsz(3)-scrsz(1) 2*scrsz(4)/3],...
            'color',[1 1 1]);
        mov = length(C(1,:));
        pltcol = [0 0 0];       
    else
%         figure('color',[1 1 1],...
%             'position',[500 50 1.5*[420 594]])       
        mov = length(C(1,:));
        pltcol = [0 0 0];
    end
else
    h = gcf;    
    pltcol = abs(get(h,'color')-[1 1 1]);   
end



if isempty(nums)
     nums = 1:length(Df);
end

%% spatial components

    
if isempty(ref) && horizontal == 1
    subplot(1,nxplot,1)
    d1 = options.d1;
    d2 = options.d2;
    if max(Cn(:))>1
        imshow(Cn/max(Cn(:))); 
    else
        imshow(Cn);
    end
    hold on;

    if mov < length(C(1,:))
        a = zeros(d1,d2);
        for i = 1:length(A(1,:))  
            aa = full(reshape(A(:,i),d1,d2));
            a = a + aa(:,:,[1 1 1])/max(max(aa))*(C(i,mov)-min(C(i,:)))/(max(C(i,:)-min(C(i,:))));
        end
    else
        a = full(reshape(sum(A(:,:),2),d1,d2));
        a = a(:,:,[1 1 1])/max(max(a));
    end
    a(:,:,[1 3]) = 0;
    fa = imshow(a);
    set(fa, 'AlphaData',a(:,:,2))

    if exist('AA','var')
    a = full(reshape(sum(AA(:,:),2),d1,d2));
    a = a(:,:,[1 1 1])/max(max(a));
    a(:,:,[2 3]) = 0;
    faa = imshow(a);
    set(faa, 'AlphaData',a(:,:,1))
    end
    hold off
elseif horizontal == 1
%     subplot(1,nxplot,1)
    axes('position',[.01 .01 .45 1])
    imshow(ref)           
end

%% Bulk Signal plot
    
if isfield(caim,'bulk')
%     subplot(6*nyplot,nxplot,nxplot-posplot)
    a1 = axes('position',[.51 .83 .45 .1]);
    if isfield(caim.bulk,'trace')
        bulk = smooth(caim.bulk.trace(:,2),10);
    else
        bulk = smooth(caim.bulk.traceMEC(1,:),10);
    end

    plot(tsscn(1:end),bulk,...
        'color',[1 0 0],...
        'linewidth',lnwd); 
    hold on
    axis tight
    axis off
    if mov < length(C(1,:))
        plot([tsscn(mov)+timeint(1) tsscn(mov)+timeint(1)],[min(bulk) max(bulk)],...
            'color',[0 1 0],...
            'linewidth',lnwd);   
    end

    set(gca,'xlim', timeint);
    hold off


    fig_pos=get(gca,'position');
    yp1=fig_pos(2);
    yp2=yp1+fig_pos(4);%/(abs(min(aa))+max(aa));%+1/(length(nums)*max(a));
    xp1=fig_pos(1); % startpoint for the scale bars

    if mov < length(C(1,:))   
    annotation('textarrow',[xp1+.095,xp1+.095],[yp2,yp2],...
       'color',pltcol,...
       'HeadStyle','none',...
       'string','MPP input Signal',...
       'fontsize',ftsz,...
       'VerticalAlignment','bottom',...  
       'HorizontalAlignment','left'); 
    end
end
%% Plot imaging Data

% subplot(6*nyplot,nxplot,nxplot*[2 3 4]-posplot)
% delete(a1.Children)
a1 = axes('position',[.51 .32 .45 .55]);

for i = 1:length(nums)

%     aa = Y(nums(i),:)/Df(nums(i));
%     aa = S(nums(i),:)/max(S(nums(i),:));

    aa = C(nums(i),:)/Df(nums(i));
    aa = aa-min(aa);
    plot(tsscn(1:end),aa(1:end)+((i-1)),...
        'color',pltcol,...
        'linewidth',lnwd);
    hold on

    if mov < length(C(1,:))
        plot([tsscn(mov)+timeint(1) tsscn(mov)+timeint(1)],[0 length(nums)],...
            'color',[0 1 0],...
            'linewidth',lnwd);   
    end
end

set(gca,'ylim', [0 length(nums)]);
set(gca,'xlim', timeint);
axis off

% scale bars

    fig_pos=get(gca,'position');
    yp1=fig_pos(2);
    yp2=fig_pos(2)+fig_pos(4)/length(nums);%/(abs(min(aa))+max(aa));%+1/(length(nums)*max(a));
    xp1=fig_pos(1)-.008;
    xp2=xp1+fig_pos(3)/length(aa)/dt*1000;


%     h1 = annotation('textarrow',[xp1,xp2],[yp1,yp1],...
%            'color',pltcol,...
%            'HeadStyle','none',...       'HeadWidth',6,...       'Headlength',12,...
%            'linewidth',lnwd);
       
h2 = annotation('textarrow',[xp1,xp1],[yp1+.16,yp1+.16],...
       'color',pltcol,...
       'HeadStyle','none',...
       'string','100% \DeltaF/F',...
       'fontsize',ftsz-2,...
       'VerticalAlignment','bottom',...  
       'TextRotation',90,...
       'HorizontalAlignment','right');

h3 = annotation('textarrow',[xp1,xp1],[yp1,yp2],...
       'color',pltcol,...
       'HeadStyle','none',...
       'HeadWidth',6,...
       'Headlength',12,...
       'linewidth',lnwd);       

% h4 = annotation('textarrow',[xp1,xp1],[yp1,yp1],...
%        'color',pltcol,...
%        'HeadStyle','none',...
%        'string','10 sec',...        
%        'fontsize',ftsz,...
%        'VerticalAlignment','top',...
%        'TextRotation',0,...
%        'HorizontalAlignment','left');

% fig_pos=get(gca,'position');
%     yp1=fig_pos(2);
    yp2=yp1+fig_pos(4);%/(abs(min(aa))+max(aa));%+1/(length(nums)*max(a));
    xp1=fig_pos(1); % startpoint for the scale bars
    
if mov < length(C(1,:))       
hcap = annotation('textarrow',[xp1+.16,xp1+.16],[yp2+.01,yp2+.01],...
       'color',pltcol,...
       'HeadStyle','none',...
       'string','Sham-Control',...
       'fontsize',ftsz,...
       'VerticalAlignment','bottom',...  
       'HorizontalAlignment','left'); 
end
hold off

 
%% Plot Belt Data

beltmov = tsscn(abs(tsscn-tsscn(mov)) == min(abs(tsscn-tsscn(mov))));
% beltend = Belt.time(abs(Belt.time-tsscn(end)) == min(abs(Belt.time-tsscn(end))));

%% Movement Plot  

% subplot(6*nyplot,nxplot,nxplot*5-posplot)    
a1 = axes('position',[.51 .12 .45 .1]);

plot(tsscn,Belt.distance/max(Belt.distance),...
        'color',[0 176/255,240/255],...
        'linewidth',lnwd);

set(gca,'ylim', [0 1]);
set(gca,'xlim', timeint);
axis off
hold on

% % Stimulation marker
% aa = find(Belt.airpuff == 1);
% for i = 1:length(aa)
%     plot([Belt.time(aa(i)) Belt.time(aa(i))],[0 1],...
%         'color',[1 0 1],...
%         'linewidth',lnwd); 
% end

% Position marker
if mov < length(C(1,:))
    plot([beltmov+timeint(1)  beltmov+timeint(1) ],[0 1],...
        'color',[0 1 0],...
        'linewidth',lnwd);   
end
hold off

% Scalebars

    
    fig_pos=get(gca,'position');
    yp1=fig_pos(2);
    yp2=yp1+fig_pos(4);%/(abs(min(aa))+max(aa));%+1/(length(nums)*max(a));
    xp1=fig_pos(1); % startpoint for the scale bars
    if mov < length(C(1,:))   
    hcap = annotation('textarrow',[xp1+.045,xp1+.045],[yp2+.01,yp2+.01],...
           'color',pltcol,...
           'HeadStyle','none',...
           'string','Position',...
           'fontsize',ftsz,...
           'VerticalAlignment','bottom',...  
           'HorizontalAlignment','left'); 
    end
    h1 = annotation('textarrow',[xp1,xp1],[yp1,yp2],...
           'color',pltcol,...
           'HeadStyle','none',...       'HeadWidth',6,...       'Headlength',12,...
           'linewidth',lnwd);
       
    h3 = annotation('textarrow',[xp1,xp1],[yp1,yp1],...
       'color',pltcol,...
       'HeadStyle','none',...
       'string','0',...
       'fontsize',ftsz-5,...
       'VerticalAlignment','middle',...  
       'HorizontalAlignment','right');    

    h4 = annotation('textarrow',[xp1,xp1],[yp2,yp2],...
       'color',pltcol,...
       'HeadStyle','none',...
       'string','1.5 m',...
       'fontsize',ftsz-5,...
       'VerticalAlignment','middle',...  
       'HorizontalAlignment','right');     
   
%% Pupil size plot

%%% a1 = axes('position',[.51 .05 .45 .1]);
%%% if isfield(Belt,'pupil')       
%%%    aa = Belt.pupil(1:length(tsscn));
%     aa(aa<.7*max(aa)) = .7*max(aa);
%%%    aa = (aa-min(aa))/max(aa-min(aa));
%%% else 
 %%%   aa = ones(1,length(Belt.time));
 %%% end 

 %%% plot(tsscn,aa,...
   %%%      'color',[.5 .5 .5],...
  %%%       'linewidth',lnwd);

 %%% set(gca,'ylim', [.7 1]);
 %%% set(gca,'xlim', timeint);
 %%% axis off
 %%% hold on

% Stimulation marker
 %%% if isfield(Belt,'airpuff')
  %%%   aa = find(Belt.airpuff.stimon == 1);
  %%%   for i = 1:length(aa)
 %%%        plot([tsscn(aa(i)) tsscn(aa(i))],[0 1],...
 %%%            'color',[1 0 1],...
  %%%           'linewidth',lnwd); 
 %%%    end
 %%% end

% Postion marker
 %%% if mov < length(C(1,:))
  %%%   plot([beltmov+timeint(1)  beltmov+timeint(1) ],[0 1],...
  %%%      'color',[0 1 0],...
  %%%      'linewidth',lnwd);   
 %%% end  
   
    
    
% Scalebars

    
 %%%fig_pos=get(gca,'position');
 %%%yp1=fig_pos(2);
 %%%yp2=yp1+fig_pos(4);%/(abs(min(aa))+max(aa));%+1/(length(nums)*max(a));
 %%%xp1=fig_pos(1); % startpoint for the scale bars
 %%%if mov < length(C(1,:))       
 %%%    hcap = annotation('textarrow',[xp1+.055,xp1+.055],[yp2+.01,yp2+.01],...
   %%%     'color',pltcol,...
  %%%      'HeadStyle','none',...
  %%%      'string','Pupil Size',...
 %%%       'fontsize',ftsz,...
 %%%       'VerticalAlignment','bottom',...  
 %%%       'HorizontalAlignment','left'); 
 %%%end

 %%%annotation('textarrow',[xp1,xp1],[yp1,yp2],...
  %%%      'color',pltcol,...
  %%%      'HeadStyle','none',...       'HeadWidth',6,...       'Headlength',12,...
  %%%      'linewidth',lnwd);
 %%%if isfield(Belt,'airpuff') 
  %%%   annotation('textarrow',[xp1+0.03,xp1+.04],[yp1-.04,yp1-.04],...
  %%%      'color',[1 0 1],...
  %%%      'HeadStyle','none',...       'HeadWidth',6,...       'Headlength',12,...
  %%%      'linewidth',lnwd,...
  %%%      'string','Airpuff',...
  %%%      'fontsize',ftsz-5,...
  %%%      'VerticalAlignment','middle',...  
  %%%      'HorizontalAlignment','left');
 %%%end

 %%%annotation('textarrow',[xp1,xp1],[yp1,yp1],...
   %%%     'color',pltcol,...
   %%%     'HeadStyle','none',...
    %%%    'string','min',...
    %%%    'fontsize',ftsz-5,...
    %%%    'VerticalAlignment','middle',...  
    %%%    'HorizontalAlignment','right');    

 %%%annotation('textarrow',[xp1,xp1],[yp2,yp2],...
  %%%      'color',pltcol,...
   %%%     'HeadStyle','none',...
  %%%      'string','max',...
  %%%      'fontsize',ftsz-5,...
  % %%%%%      'VerticalAlignment','middle',...  
 %      'HorizontalAlignment','right');    
   
%% timescalebar

lbar = 30; %length of scalebar in seconds   
% % fulllength
% tmscl = lbar*1000*(length(Belt.time)/Belt.time(end))*(fig_pos(3)/(length(Belt.time))) % legnth of the scale bar normalized to length of plot and recording
% interval
tmscl = lbar*1000*(fig_pos(3)/(timeint(2)-timeint(1)));
xp1=fig_pos(1)+fig_pos(3)-tmscl; 
xp2=xp1+tmscl;                  
  
annotation('textarrow',[xp1,xp2],[yp1,yp1],...
       'color',pltcol,...
       'HeadStyle','none',...
       'HeadWidth',6,...
       'Headlength',12,...
       'linewidth',lnwd); 
   
annotation('textarrow',[xp2,xp2],[yp1,yp1],...
       'color',pltcol,...
       'HeadStyle','none',...
       'string','30 sec',...        
       'fontsize',ftsz,...
       'VerticalAlignment','top',...
       'TextRotation',0,...
       'HorizontalAlignment','left');
    
hold off

%%
if save
    [filename, pathname] = uiputfile('*.png','Choose figure name','/media/2Photon/Martin Pofahl/');

    if ~filename == 0
        exportfig([pathname filename(1:end-4)])
    end
end
   

end
   