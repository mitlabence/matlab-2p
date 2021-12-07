function FFMPEGConvert(input_videos_folder, input_extension, output_path, output_extension)
%FFMPEGCONVERT Given a folder where videos (for example, HIK mp4 videos)
% of a given extension are found, this function converts them into a given
% extension and puts them in a folder supplied as an argument.
% 
% input_videos_folder:  path name ending with "\", where the videos are
% input_extension:  extension of videos to convert with ".", one example:
%                   ".mp4"
% output_path: path name ending with "\", where to put the converted vids.
% output_extension: into which format to convert the videos. Example:
%                   ".avi"

%% Get list of videos
vids = dir(strcat(input_videos_folder, "*", input_extension));

%% Use FFMPEG to convert proprietary mp4 to avi
for vid = vids'
    [~, vid_name, vid_extension] = fileparts(vid.name); % first is path, third is extension (".mp4"), not needed
    input_file = strcat(input_videos_folder, vid_name, vid_extension); % not sure if fileparts() just above detects the path! Hence the use of VIDEOS_FOLDER
    output_file = strcat(output_path, vid_name, output_extension); 
    system(sprintf('ffmpeg -err_detect ignore_err -i %s -c copy %s', input_file, output_file));
end

%% Display text on completion

%converted_vids = dir(strcat(output_path, "*", output_extension)); % list
%of all converted videos. Returning it could be optional, but not necessary

disp(fprintf("FFMPEG conversion from %s to %s complete", input_extension, output_extension));

end

