function frame = GetFrameForPreview(video_fpath,frame_number, temp_path)
%GETFRAME From the video at video_fpath, get the given frame frame_number (indexing starts with 0) from it. Temporarily
%save image in temp_path, delete at the end of this function.
%   video_fpath: Full path of video to get a preview frame from.
%   frame_number: The index of the frame to be grabbed (0 is the first
%   frame)
%   temp_path: Full path of temporary folder where an image can be saved
%   and then deleted during the execution of this funcition. May or may not
%   end with "\" (i.e. "C:\Images" is valid and equivalent to the also
%   correct "C:\Images\".

%TODO: error handling if arguments are missing or invalid

%% Create output file path
output_fpath = strcat('preview_', GetTimestampForFilename(), '.png');  %put image in temp folder, delete automatically later
if(~endsWith(temp_path, '\'))
    output_fpath = strcat('\', output_fpath);
end
output_fpath = strcat(temp_path, output_fpath);

%% Get frame
system(sprintf('ffmpeg -i %s -vf "select=eq(n\\, %i)" -vframes 1 %s', video_fpath, frame_number, output_fpath));
frame = imread(output_fpath);

%% Clean up (delete temporarily saved image)
if exist(output_fpath, 'file')==2
    delete(output_fpath);
end
end

