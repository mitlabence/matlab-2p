%Dependencies:
% ffmpeg: need to be added to Path. Used as command-line command (ffmpeg
% ... in powershell)
% videoMultiCrop:   A modified version is used by this code. The original:
%                   https://de.mathworks.com/matlabcentral/fileexchange/34845-crop-video-to-sub-elements
%                   For selecting parts of the video (mouse cages), to get
%                   better signal/noise ratio and able to observe
%                   individual mice from the movement perspective.
% Signal Processing Toolbox
% Possibly: Computer Vision Toolbox (not sure if it is used)

%% Testing part, only for developing uses
%addpath 'D:\PhD\Matlab Scripts\extern\videoMultiCrop';

%input_fpath = 'D:\PhD\testvid\';
%input_fname = 'hiv00008.mp4';
%output_fpath = input_fpath;
%output_fname = 'converted_vid.avi';

%input_file = strcat(input_fpath, input_fname);
%output_file = strcat(output_fpath, output_fname);

%% Add external libraries
%addpath 'D:\PhD\Matlab Scripts\extern\videoMultiCrop'; % add
%videoMultiCrop. Not needed because I use modified videoMultiCrop

%% Select the files

VIDEOS_FOLDER = 'D:\PhD\Useful_Video_Monitoring\';  % Select folder where videos (.mp4, HIKVision h.264+ codec) are located
                                                    % with \ as last character
OUTPUT_PATH = 'D:\PhD\Useful_Video_Monitoring\output\';  % store converted avi videos
TEMP_PATH = 'D:\PhD\Useful_Video_Monitoring\temp\';  % store individual mouse videos
RESULTS_PATH = 'D:\PhD\Useful_Video_Monitoring\results\';

INPUT_EXTENSION = '.mp4';
OUTPUT_EXTENSION = '.mp4';

MOUSE_NAMES = ["T532", "T349"]; %if we used '', it would be concatenated...

%TODO: create folders above if they do not exist

%% Convert with FFMPEG (not needed anymore!)
%FFMPEGConvert(VIDEOS_FOLDER, INPUT_EXTENSION, OUTPUT_PATH, OUTPUT_EXTENSION);

%% Get rectangles (ROI) for mice
%TODO: ffmpeg -i in.avi -filter:v "crop=out_w:out_h:x:y" out.avi, see https://video.stackexchange.com/questions/4563/how-can-i-crop-a-video-with-ffmpeg
%much faster than what I currently have! Also no need to convert to avi, if
%I can set the rectangles in matlab on the original videos!
%Able to check the rectangle using: ffplay -i input.mp4 -vf "crop=in_w:in_h:x:y"

videosList = dir(strcat(VIDEOS_FOLDER, '*.mp4'));
%use the first video found to set the rectangles to crop
%TODO: save image with rectangles! But only after bug with last mouse name
%not showing up has been solved
rects = GetCropRectanglesFromPreview(strcat(videosList(1).folder, "\", videosList(1).name), 0, MOUSE_NAMES, TEMP_PATH);
save(strcat(RESULTS_PATH, 'rects_', GetTimestampForFilename(), '.mat' ), 'rects'); % save rects variable for later use
disp(rects);

nSubVideos = size(rects); % number of actual rectangles selected.
nSubVideos = nSubVideos(2); % rects is 1 x n cell, we need n

%load('D:\PhD\testvid\test\results\rects.mat'); % load rects variable
%TODO: find a way to choose whether to load rects variable or create new
%ones (when perspective changes)

%% Get optical flow time series and process it
for vid = videosList'
    outputVideos = cell(1, nSubVideos);
    inputVidFullPath = strcat(vid.folder, "\", vid.name);
    splitVidNameExtension = split(vid.name, "."); %get original video name without extension (1st element of this split)
    %out_vids = videoMultiCrop('inVideo', strcat(vid.folder, '\', vid.name), 'rectPosCell', rects, 'outputFolder', TEMP_PATH); %TODO: change output videos output folder
    %Process in two loops: first, create cropped videos and save them.
    %Second, create the flow time series.
    for iSubVideo = 1:nSubVideos %
        outputVidFullPath = strcat(OUTPUT_PATH, splitVidNameExtension{1}, "_", MOUSE_NAMES{iSubVideo}, "_", GetTimestampForFilename(), OUTPUT_EXTENSION );
        CropVideoToRectangle(inputVidFullPath, rects{iSubVideo}, outputVidFullPath);
        outputVideos{iSubVideo} = outputVidFullPath;
    end
    for iSubVideo = 1:nSubVideos
        subTitle = strcat(splitVidNameExtension(1), " ", MOUSE_NAMES(iSubVideo)); %This should be used when saving figure and csv
        disp("Now processing %s", MOUSE_NAMES{iSubVideo} ); %FIXME: %s is displayed, not the number. sprintf needed?
        croppedVideo = outputVideos{iSubVideo};
        [t, flow] = GetOpticalFlow(croppedVideo); %TODO: more informative text, not only 
        flow_filtered = RemoveArtificialSpikes(flow);
        h = figure;
        plot(t, flow_filtered);
        ylim([-max(flow_filtered), 5*max(flow_filtered)]);
        title(subTitle);
        %TODO: check naming convention, and title (should contain video
        %name)
        savefig(h, strcat(RESULTS_PATH, splitVidNameExtension(1), '_', MOUSE_NAMES(iSubVideo), "_", GetTimestampForFilename(), '.fig'));
        writematrix([t; flow'; flow_filtered'], strcat(RESULTS_PATH, splitVidNameExtension(1), '_', MOUSE_NAMES(iSubVideo), "_", GetTimestampForFilename(), '.csv')); %flow and flow_filtered need to be transposed
        close(h);
    end
end
%NEXT: try to compare ANfall with MagicalFilter-filtered time series.
%Probably either object tracking (try to find the one I heard about in one
%of the tuesday seminars!) or ROI-separation (check the upper part of the
%cage for jumping mice to find at least those kinds of seizures).

%TODO: VERY IMPORTANT! Make sure that both the non-filtered and filtered
%flow are in the csv saved (i.e. flow and flow_filtered are correct), and
%that the flow is not changed in RemoveArtificialSpikes().

%TODO: length(y) gives a single number, no need to use two-line size() then
%(2) or (1) if row or column! i.e. if size(y) = 1x300, then length(y) =
%300. No need for l = length(y), l = l(2)!

%TODO: all movement is visible as signal. How about checking the upper
%half? Only when the mouse has a jumpy seizure would we expect to see
%changes there.

%TODO: somehow get where the iframes are in the video? And replace those by
%the averages of the neighbors. Maybe that could be an effective filtering.

%TODO: find a way to make crop selection data accessible later on, like
%give names of mice first, then run crop function, where 1st rectangle will
%be 1st mouse etc. Then at the end, show each signal in same plot (hold
%on), as well as in a different figure the first frame of the video with
%the rectangles with corresponding color and name. Save the image and the
%signals.


