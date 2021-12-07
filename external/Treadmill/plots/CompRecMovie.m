function CompRecMovie(caim,scn,nam)
%% Network movie
timeint = [16 20]; %in min
timeint = timeint*60*1000;

options.sframe = find(abs(scn.tsscn-timeint(1))==min(abs(scn.tsscn-timeint(1))),1);						% user input: first frame to read (optional, default 1)
options.num2read = find(abs(scn.tsscn-timeint(2))==min(abs(scn.tsscn-timeint(2))),1)-options.sframe+1;					% user input: how many frames to read   (optional, default until the end)
options.crop = [30 30 30 30];
options.numchan  = 2;

[Data,Datared] = readdata(nam,options);


%% Average movie
avfc = 5;
movlg = floor(size(Data,3)/avfc);
Datamean = zeros(size(Data,1),size(Data,2),movlg);
Dataredmean = zeros(size(Data,1),size(Data,2),movlg);
for i = 1:movlg
    Datamean(:,:,i) = mean(Data(:,:,(i-1)*avfc+1:i*avfc),3);
    Dataredmean(:,:,i) = mean(Datared(:,:,(i-1)*avfc+1:i*avfc),3);
end 
Datamean = mat2gray(Datamean);
Dataredmean = mat2gray(Dataredmean);
%% adjust contrast
contrastred = stretchlim(Dataredmean(:,:,1));
contrastred(1)  = contrastred(1)/2;
contrastred(2)  = contrastred(2)*1.1;
contrastgreen = stretchlim(Datamean(:,:,1));
contrastgreen(1)  = contrastgreen(1)/2;
contrastgreen(2) = contrastgreen(2)*1.3;
%% initiate video
vidObj = VideoWriter([nam(1:end-4) ' recon.avi']);
set(vidObj,'FrameRate',15);
open(vidObj);
%% set figure (with color snd size)
h = figure('color',[0 0 0],...
    'Position',[300 250 1400 570],...
    'visible','off');

%%
j=1;
for i = 1:movlg
    frame = imadjust(Dataredmean(:,:,i),contrastred);
    frame(:,:,2) = imadjust(Datamean(:,:,i),contrastgreen);
    frame(:,:,3) = 0; 
    %%
    h = makeframe(caim,frame,i+floor(options.sframe/avfc),avfc);
    currFrame = getframe(h);
    writeVideo(vidObj,currFrame);
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

function h = makeframe(caim,frame,mov,avfc)
%%

h = gcf; 
options = caim.options;
d1 = options.d1;
d2 = options.d2;
A = caim.A;
A = reshape(full(A),d1,d2,size(A,2));
% A = mat2gray(A);
C = caim.C;

if isempty(frame)
    Cn = caim.Cn;
    Cn = cat(3, Cn, Cn, Cn);
else
    Cn = frame;
end
% if max(Cn(:))>1
%    Cn = (Cn-min(Cn(:)))./(max(Cn(:)-min(Cn(:))));
% end
graymode = 0; %set this to 1 if you dont want colored components
if graymode == 1
    AA = zeros(d1,d2);
    for j = 1:size(C,1)
        Atemp = zeros(d1,d2);
        for i = 1:5
            time = mov*avfc+i-1;
            Atemp = Atemp+A(:,:,j)*C(j,time);        
    %         Atemp = Atemp/30;
        end
        AA = AA+Atemp/400;
    end
    % AA = mat2gray(AA);
    AA = cat(3,AA,AA,AA);
else
    AA = zeros(d1,d2,3);
    % you can add more colors here if you want
    compcolor = [1 0 0;1 0 1; 1 1 0;0 1 0;0 0 .8;.5 0 1;1 .5 0;0 1 1;1 0 .5;0 .5 0;.8 .5 0;0 .5 .5;.5 1 .5;];
    j = 1;
    for k = 1:size(C,1)
        Atemp = zeros(d1,d2);
        for i = 1:5
            time = mov*avfc+i-1;
            Atemp = Atemp+A(:,:,k)*C(k,time);        
    %         Atemp = Atemp/30;
        end
        Atemp = cat(3,Atemp,Atemp,Atemp);
        for jj = 1:3
            Atemp(:,:,jj) = Atemp(:,:,jj).*compcolor(j,jj);    
        end
        AA = AA+Atemp/300;
        j = j+1;
        if j == size(compcolor,1)
            j = 1;
        end
    end
end

% concatenate picture and components
c = [Cn AA];
imshow(c);

end