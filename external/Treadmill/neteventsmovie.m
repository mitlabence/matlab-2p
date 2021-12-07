function neteventsmovie(caim,nam)

%% Network movie

 nam = '/media/2Photon/Nicola/339-342-349/130919/339.130919/zone1/M339.130919.1538 out.tif';
options.sframe=1;						% user input: first frame to read (optional, default 1)
options.num2read=[];					% user input: how many frames to read   (optional, default until the end)
options.crop = [20 20 20 20];
options.numchan  = 2;

[Data,~] = readdata(nam,options);
%%
avfc = 5;
Z = zeros(size(Data,1),size(Data,2),floor(size(Data,3)/avfc));
for i = 1:size(Z,3)
    Z(:,:,i) = mean(Data(:,:,(i-1)*avfc+1:i*avfc),3);
end 

Z = Z - min(Z(:)); 
% Cn = (Cn-min(Cn(:)))./(max(Cn(:)-min(Cn(:))));
Z = Z/max(max(mean(Z,3)));
movlg = floor(size(Z,3));
%%
vidObj = VideoWriter([nam(1:end-4) '.avi']);
set(vidObj,'FrameRate',15);
open(vidObj);
%% 
scrsz = get(0,'ScreenSize');
h = figure('color',[0 0 0]);
% 'Position',[scrsz(1) scrsz(4)/4 scrsz(3)-scrsz(1) 2*scrsz(4)/3]);

set(h,'visible','off');  

j=0;

%%
for i = 1:movlg
    
    h = makeframe(caim,Z(:,:,i),i);
    currFrame = getframe(h);
    writeVideo(vidObj,currFrame);
    pause(.001)
    if j == 50
        disp([num2str(i) '/' num2str(movlg)])
        j = 0;
    end
    delete(allchild(h))
    j = j+1;
    
end
delete(h)
close(vidObj)
disp('fertig')




end

function h = makeframe(caim,frame,mov)
A = caim.A;
C = caim.C;
Df = caim.Df;
options = caim.options;
d1 = options.d1;
d2 = options.d2;
netID = caim.network.netID;
netpos = caim.network.netpos;

if isempty(frame)
    Cn = caim.Cn;
else
    Cn = frame;
end
%%

% if max(Cn(:))>1
%    Cn = (Cn-min(Cn(:)))./(max(Cn(:)-min(Cn(:))));
% end

Cn = cat(3, Cn, Cn, Cn);

  %%  
% if isempty(mov)
%     scrsz = get(0,'ScreenSize');
%     h = figure('Position',[scrsz(1) scrsz(4)/4 scrsz(3)-scrsz(1) 2*scrsz(4)/3],...
%         'color',[1 1 1]);
%     mov = length(C(1,:));
%     pltcol = [0 0 0];
% else
    h = gcf;
    set(h,'color',[0 0 0]);
    pltcol = [1 1 1];
% end

%% plot networkevents-components

b = full(reshape(sum(A(:,netID{1}),2),d1,d2));
b = b(:,:,[1 1 1])/max(max(b));
b(:,:,:) = 0;

j = 0;
time = mov * 5;
for i = 1:length(netID)
    if netpos(i) < time    
        a = full(reshape(sum(A(:,netID{i}),2),d1,d2));
        a = a(:,:,[1 1 1])/max(max(a));
        if j == 0
            a(:,:,[1 3]) = 0;
            j = j+1;
        elseif j == 1
            a(:,:,[1 2]) = 0;
            j = j+1;
        elseif j == 2
            a(:,:,[2 3]) = 0;
            j = j+1;
        elseif j == 3
            a(:,:,[1 ]) = 0;
            j = j+1;
        elseif j == 4
            a(:,:,[2 ]) = 0;
            j = j+1;
        elseif j == 5
            a(:,:,[3 ]) = 0;
            j = 0;
        end  
        b = b + a;  
    end  
%     imshow(b);
%     pause(.5)
end

c = [Cn b];
faa = imshow(c);
end