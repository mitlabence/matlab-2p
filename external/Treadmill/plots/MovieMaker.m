%% The Presentation Movie Maker
% Choose the path to the data file you want to show
pathname = '/media/2Photon/Nicola/VIDEOpaper/349ica/202020/';

trial = 1;
[files,~,prefiles] = FindDataSets(pathname);
disp(['load ' prefiles(trial).name])
load(prefiles(trial).name)

%% This writes a data movie 
% make sure to set the path to the right tif
MovieComp(caim,scn,[pathname files(1).name(1:end-6) '.tif'])
%% This writes a component reconstruction movie
CompRecMovie(caim,scn,[pathname files(1).name(1:end-6) '.tif'])