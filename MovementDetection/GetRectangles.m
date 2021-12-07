function [rectPosCell] = GetRectangles(varargin)
%GETRECTANGLES Summary of this function goes here
%   Detailed explanation goes here
%% Default params
inVideo = [];
iFrames = [];
selectFrame = [];
rectPosCell = {};
miceList = [];

%% Load uses params, overifding default ones
% assert(nargin > 0, 'Inputs missing');
% automatically get all input pairs (name, value) and store in local vairables
for iArg=1:2:length(varargin)
    eval([varargin{iArg},'=varargin{iArg+1};']);  % TODO get read of EVAL.
    %    assignin_value(varargin{iArg}, varargin{iArg+1});
end

%% Prepare input video object
videoFormats= VideoReader.getFileFormats();
videosFilesExtList={videoFormats.Extension};
inVideo=filesFullName(inVideo, videosFilesExtList);

[videoPath, videoName, videoExt] = fileparts(inVideo);
ObjInVideo = VideoReader([videoPath, filesep, videoName, videoExt]);

% in case videoPath was missing- get it from ObjInVideo fields
inVideo=strcat(ObjInVideo.Path, filesep, ObjInVideo.Name);
[videoPath, videoName, videoExt] = fileparts(inVideo);
outVideoExt=videoExt;

nFrames=ObjInVideo.NumberOfFrames;
if exist('iFrames', 'var')~=1 || isempty(iFrames)
    iFrames = 1:nFrames;
else
    isLegalFrame = iFrames>0 & iFrames<=nFrames;
    if any(~isLegalFrame)
        iFrames = iFrames(isLegalFrame);
    end
end
if exist('selectFrame', 'var')~=1 || isempty(selectFrame)
    selectFrame = iFrames(1);
end

n_rect = 0;
N_mice = size(miceList);

%% Set the sub video coordinates by dragging a rectangle over a video frame
if isempty(rectPosCell)
    iSubMovie=0;
    anotherRect='One more';
    
    prompt = {'Enter frame index.'};
    dlgTitle = 'Choose frame';
    def = {sprintf('%d', selectFrame)};
    
    figureH=figure;
    anotherFrame='Another frame';
    
    firstTime=true;
    while strcmpi(anotherFrame, 'Another frame') 
        if firstTime
            firstTime=false;
        else
            selectFrame = str2double( inputdlg(prompt, dlgTitle, 1, def) );
        end
        
        currFrame=read(ObjInVideo, selectFrame);
        imshow(currFrame);
        title(sprintf('Movie frame N'' %d', selectFrame), 'FontSize', 14,...
            'Color', [1, 0, 0]);
        anotherFrame = questdlg({'Does this frame siuts your needs? ',...
                'Press ''Another frame'', to select another Frame.',...
                'Press ''Ok'' to use this frame.'},...
                'Frame selection',...
                'Another frame', 'Ok','Ok');
    end
    
    textH=text(  size(currFrame,2)/2 ,size(currFrame,1)+30,...
        sprintf('Choose sub movie area for mouse %s. \nDouble click when done.', miceList(n_rect+1)), 'FontSize', 14,...
        'Color', [0,0,1], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    
    while strcmpi(anotherRect,'One more')
        iSubMovie=iSubMovie+1;
        delete(textH);
        textH=text(  size(currFrame,2)/2 ,size(currFrame,1)+30,...
            sprintf('Choose sub movie area for mouse %s. \nDouble click when done.', miceList(n_rect+1)), 'FontSize', 14,...
            'Color', [0,0,1], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
        
        rectH = imrect;
        textROIPosChanged(getPosition(rectH));
        id = addNewPositionCallback(rectH, @textROIPosChanged); % display selection info
        rectPosCell{iSubMovie}=wait(rectH); %getPosition(rectH); % [xmin ymin width height]
        delete(rectH);
        
        % verify legal values- integers which fit currFrame dimentions
        rectPosCell{iSubMovie}=round(rectPosCell{iSubMovie});
        rectPosCell{iSubMovie}(1:2)=max(rectPosCell{iSubMovie}(1:2),1);
        
        rectPosCell{iSubMovie}(3:4)=min(rectPosCell{iSubMovie}(3:4),...
            [size(currFrame,2)-rectPosCell{iSubMovie}(1),...
            size(currFrame,1)-rectPosCell{iSubMovie}(2)]);
        
        rectangle( 'Position', rectPosCell{iSubMovie}, 'Curvature', [0, 0],...
            'EdgeColor', rand(1, 3), 'LineWidth', 2 );
        
        n_rect = n_rect + 1;
        
        if (~isempty(miceList) && n_rect >= N_mice(2))
            break
        end
        anotherRect = questdlg({'Care to mark another sub video rectangle? ',...
            'Press ''One more'', to mark another sub video.',...
            'Press ''Finish'' to finish choosing input files.'},...
            'Outputs files area selection',...
            'One more', 'Finish','Finish');
    end
    delete(textH);
    close(figureH);
end % if isempty(rectPosCell)
if ~iscell(rectPosCell)
    rectPosCell={rectPosCell};
end

%% service subfunction
function textROIPosChanged(pos)
pos=round(pos);
poseMsg=sprintf('Rectangle dimensions: [left, bottom, heigh, width]=[%d, %d, %d, %d]',...
    pos(1), pos(2), pos(3), pos(4));

title(poseMsg, 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [0.7, 0.9, 0.7]);
