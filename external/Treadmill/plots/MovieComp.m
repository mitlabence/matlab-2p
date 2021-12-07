function []= MovieComp(caim,scn,pathname)
%% load data
addpath(genpath('/media/2Photon/Matlab code/Nicos Source Code/Treadmill/plots/'));
addpath(('/media/2Photon/Matlab code/Nicos Source Code/'));
% pathname = ('/media/2Photon/Nicola/VIDEOpaper/339.13/');
%%
timeint = [16 20]; %in min
timeint = timeint*60*1000;

options.crop = [0 0 0 0];
options.sframe = find(abs(scn.tsscn-timeint(1))==min(abs(scn.tsscn-timeint(1))),1);						% user input: first frame to read (optional, default 1)
options.num2read = find(abs(scn.tsscn-timeint(2))==min(abs(scn.tsscn-timeint(2))),1)-options.sframe+1;					% user input: how many frames to read   (optional, default until the end)
options.numchan  = 2;
%%
[Data,Datared] = readdata(pathname,options);

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
%%
vidObj = VideoWriter([pathname(1:end-8) ' data.avi']);
set(vidObj,'FrameRate',15,...
    'Quality',30);
open(vidObj);
%%
h = figure('Position',[300 200 1400 570],...
    'color',[1 1 1],...
    'renderer','painters',...
    'visible','off');

%%
contrastred = stretchlim(Dataredmean(:,:,1));
contrastred(1)  = contrastred(1)/2;
contrastred(2)  = contrastred(2)*1.1;
contrastgreen = stretchlim(Datamean(:,:,1));
contrastgreen(1)  = contrastgreen(1)/1;
contrastgreen(2) = contrastgreen(2)*2.2;
%%
j = 1;
for i = 1:movlg
    %%    
    ref = imadjust(Dataredmean(:,:,i),contrastred);
    ref(:,:,2) = imadjust(Datamean(:,:,i),contrastgreen);
    ref(:,:,3) = 0; 
    ref(:,:,1) = 0; 
%     imshow(ref)
    %%
    PlotComponentsBehave(caim,1:35,[],i*avfc,scn,ref,timeint,1)
    %%
    currFrame = getframe(h);
    writeVideo(vidObj,currFrame);
    pause(.01)
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
