function outVideos=videoMultiCrop(varargin)
%% videoMultiCrop
% This function divides the input video to several parts/sub videos.
%
%% Syntax
%  outVideo=videoMultiCrop('inVideo',videoFile,'hAppliedFunc',@function,...
%     'flagConcat',true,'compression','NONE','PARAMS',param1,param2...);
%  outVideo=videoMultiCrop('inVideo',videoFile,'hAppliedFunc',@function,...
%     'PARAMS',param1,param2...);
%
%% Description
% This functions allows the user to crop the input movie to several (unlimited) fragments
% of his choice. The fragment dimensions are user defined either via functions input or
% via GUI control (using the imrect function). The function can be thought of as the
% complementary of my concatVideo2D, (although this isn't quiet accurate). Same results
% can be reached via my apply2VideoFrames.m function
% (http://www.mathworks.com/matlabcentral/fileexchange/32351-apply2videoframes) when
% Matlab imcrop function is utilized, but when you need to get several sub videos at a
% time- this is easier and faster. Another way to reach this functionality in a single
% line command would be through the combination of apply2VideoFrames.m with funXapply.m
% (http://www.mathworks.com/matlabcentral/fileexchange/33778-apply-multiple-functions-in-a-single-function)
% but that would be a bit complicated.
%
%% Input arguments (defaults exist):
%   inVideo- input video file name
%   iFrame- an array of integers- a list of rames to be used.
%   rectPosCell-  cell array with coordinates [left, bottom, heigth, width] of each sub
%   	video.
%   compression-  compression method used. 
%   outVideos-    cell array including the the names of the final video
%     files(path+name+extension). Cell array length is assumed to fit length of
%     rectPosCell array. Note that you can choose a file type different from the input
%     file extension.
%   selectFrame-  frame which will be presented for ROI selection.
%   outputFolder - output path: 'C:\Videos' is valid, 'C:\Videos\' is also
%   valid. BUT: do NOT use ""! That is string array, not string
%
%% Output arguments
%   outVideoFile-    a cell array with the names if newly created video files
%     (path+name+extension).
%
%% Issues & Comments
% -The duality of avifile and VideoWriter makes the code cumbersome. Unfortunately no
%     single function capable of encoding all files, with all extensions  is available.
%     Another issue is that VideoWriter and VideoReader/mmreader mean different things
%     under the VideoFormat concept. IMHO until those are solved, use avifile with .AVI
%     files...
%
% -.mpg file extentin writing results in error (can be fixed if all outputs videos are
%     AVI)
%
%% Example
% matlabVfile='rhinos.avi'; % 'xylophone.mpg'; % no other input video files in Matlab default path was found
% outVideoFiles=videoMultiCrop('inVideo', matlabVfile, 'rectPosCell', {[10,10,120,180], [30,40,70,100]});
% implay(outVideoFiles{1});
% implay(outVideoFiles{2});
% implay(matlabVfile);
%
%
%% See also
%  - concatVideo2D      % http://www.mathworks.com/matlabcentral/fileexchange/33951-concatenate-video-files-subplot-style
%  - apply2VideoFrames  % http://www.mathworks.com/matlabcentral/fileexchange/32351-apply2videoframes
%
%% Revision history
% First version: Nikolay S. 2010-04-12.
% Last update:   Nikolay S. 2012-05-15.
%
% *List of Changes:*
% 2012-05-15
% - Video file is chosen via Files Explorer, if none was specified via function input
%     arguments.
% 2012-03-14
% - Some minor bugs & documentation fixed.
% - Frame chosie dialog added for GUI option of RECT choise.
%
% 2011-11-28 (inheritage from apply2VideoFrames):
% - VideoWriter added to replace avifile for non AVI files.
% - VideoReader replaced mmreader.
%


%% Default params
inVideo = [];
iFrames = [];
selectFrame = [];
outVideos = {};
rectPosCell = {};
outputFolder = '';

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

% Good time to set output folder!
if strcmp(outputFolder, '')
    outputFolder = videoPath;
    disp('No output folder given. Using input folder.');
end
% Check if output Folder ends with "\", and delete it if yes
if endsWith(outputFolder, '\')
    outputFolder = outputFolder(1:end-1);
end 
fprintf(1, 'Output folder for cropped videos is %s\n', outputFolder); %should display the output folder correctly

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


%% Set the sub video coordinates by draggin a rectangle over a video frame
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
        sprintf('Choose sub movie area. \nDouble click when done.'), 'FontSize', 14,...
        'Color', [0,0,1], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    
    while strcmpi(anotherRect,'One more')
        iSubMovie=iSubMovie+1;
        
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
        anotherRect = questdlg({'Care to mark another sub video rectangle? ',...
            'Press ''One more'', to mark another sub video.',...
            'Press ''Finish'' to finish choosing input files.'},...
            'Outputs files area selection',...
            'One more', 'Finish','Finish');
    end
    delete(textH);
    % close(figureH); uncomment this if want to close preview once all
    % rectangles selected
end
if ~iscell(rectPosCell)
    rectPosCell={rectPosCell};
end
nSubVideos=length(rectPosCell); % total number of coordinates vectors

if isempty(outVideos) % generate outVideos name if user didn't supply one
    outVideos=cell(nSubVideos, 1);
end

ObjOutVideo(nSubVideos, 1) = struct('VideoWriter', [], 'iCols', [], 'iRows', []);
for iSubMovie=1:nSubVideos
    % for each movie- if name is mising generate one automatially
    if isempty( outVideos{iSubMovie} )
        outVideos{iSubMovie}=[outputFolder, filesep, 'subV_', num2str(iSubMovie),...
            videoName, outVideoExt]; %[videoPath, ... originally
    end


    if exist('compression','var') == 1 % Use user specified compression
        ObjOutVideo(iSubMovie).VideoWriter = VideoWriter( outVideos{iSubMovie}, compression );
    else
        ObjOutVideo(iSubMovie).VideoWriter = VideoWriter( outVideos{iSubMovie} );
    end
    open( ObjOutVideo(iSubMovie).VideoWriter );
    
    posVec=rectPosCell{iSubMovie};
    
    iRows = posVec(2)+(1:posVec(4));
    isLegal = iRows>0 & iRows<=ObjInVideo.Height;
    if any(~isLegal)
        iRows=iRows(isLegal);
    end
    ObjOutVideo(iSubMovie).iRows = iRows;
    iCols = posVec(1)+(1:posVec(3));
    isLegal = iCols>0 & iCols<=ObjInVideo.Width;
    if any(~isLegal)
        iCols=iCols(isLegal);
    end
    ObjOutVideo(iSubMovie).iCols = iCols;

end % for iSubMovie=1:nSubVideos

fprintf('\nVideo %s slicing to parts started: %s \n', [videoName, videoExt],...
    datestr(now,'dd-mmm-yyyy HH:MM:SS'));


nFramesList = length(iFrames);
for iFrame = 1:nFramesList % loop through all frames
    currFrame=read( ObjInVideo, iFrames(iFrame) );
    for iSubMovie=1:nSubVideos % Crop all sub videos from each frame and store them.
        % This way we read input movie only once. Assuming written movies volume sum
        % is smaller the read movie volume, this ordering should be slightly faster.
        % in case of overlaping regions, changing loops orders may be better option.
        
        subFrame=currFrame( ObjOutVideo(iSubMovie).iRows, ObjOutVideo(iSubMovie).iCols, : );

        writeVideo(ObjOutVideo(iSubMovie).VideoWriter, subFrame);
    end % for iSubMovie=1:nSubVideos % cut out all sub videos from each frame, and save
    
    waitbarTimeRemaining2(iFrame/nFramesList, 'Cropping videos');
end % for frame_ind = 1 : nFrames % loop through all frames

for iSubMovie=1:nSubVideos
    close( ObjOutVideo(iSubMovie).VideoWriter );
end


%% service subfunction
function textROIPosChanged(pos)
pos=round(pos);
poseMsg=sprintf('Rectangle dimentions: [left, bottom, heigh, width]=[%d, %d, %d, %d]',...
    pos(1), pos(2), pos(3), pos(4));

title(poseMsg, 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [0.7, 0.9, 0.7]);