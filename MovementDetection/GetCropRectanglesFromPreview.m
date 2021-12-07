function rectCoords =  GetCropRectanglesFromPreview(video_fpath,frame_number, labels_list, temp_path)
%GETCROPRECTANGLESFROMPREVIEW Allows the user to open a frame (frame_number
%index) from the video at video_fpath, and define as many rectangles as
%labels in labels_list. The frame is temporarily saved (and deleted before 
%the function returns) to the directory temp_list
%   
%   video_fpath: full path of input video file
%   frame_number: the index of the frame in the video to open (starts with
%   0)
%   labels_list: form of ["a", "b", "c", "d"], i.e. 1 x n list of strings.
%   Use " instead of '!
%   temp_path: Full path (can end with or without "\") to folder where the
%   frame is temporarily saved as a .png image.


%% Find a good frame to draw the rectangles on
frame = GetFrameForPreview(video_fpath,frame_number, temp_path);
rectCoords = {};
len = size(labels_list);
len = len(2); % labels_list is 1 x n list of strings

%TODO: use vararg, inputVideo, frameNumber, tempPath, numberOfRectangles
%TODO: Use questdlg to get to a good frame

%% Draw the rectangles
fi = figure;
imshow(frame);
for rectangle_index = 1:len
    roi = drawrectangle; %drawrectangle
    addlistener(roi, 'ROIMoved', @showpos);
    wait(roi); %wait until user double-clicks on rectangle, making the selection final
    %save coordinates before deleting drawrectangle object
    coords = roi.Position;
    disp(coords);
    %make roi unmovable here, i.e. a constant rect.
    %delete(roi);
    roi.InteractionsAllowed = 'none';
    roi.Label = labels_list(rectangle_index);
    rectCoords{rectangle_index} = round(coords); % want to have integer coordinates -> round
    %FIXME: last rectangle does not have label
end
%% Return the rectangle list
%close(fi);
end

function showpos(src, evt)
    evtname = evt.EventName;
    switch(evtname)
        case{'ROIMoved'}
            title(mat2str(evt.CurrentPosition()), 'FontSize', 12, 'FontWeight', 'bold', 'BackgroundColor', [0.7, 0.9, 0.7]);
 
    end
end