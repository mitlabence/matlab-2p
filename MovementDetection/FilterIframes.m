function filteredTimeSeries = FilterIframes(timeSeries, iFramesInterval)
%FILTERIFRAMES Summary of this function goes here
%   timeSeries: 1 x n array of double values
%   Detailed explanation goes here
%   iFramesInterval is 50 in the HikVision recordings.
%   Checking I-Frame interval is not trivial. Useful tools are ffmpeg and
%   ffprobe. I think both are installed when ffmpeg is installed.
%   One can use 
%   ffprobe -show_frames -select_streams v:0 input_vid.mp4 > ./frames_info.txt
%   to output data about each frame into frames_info.txt in the current
%   folder. This takes a while and results in a huge txt file. It is enough
%   to wait a few seconds, then use Ctrl+C to stop the process. The txt
%   file should be still small, not complete, but contain a few
%   ten-thousand frames which is more than enough to check the period of
%   the I-Frames. These frames have the property pict_type=I, as opposed to
%   B or P. These frames are explicitely saved, not forward-predicted (like
%   P) or bi-directionally predicted (B), hence they are the outliers.

len = size(timeSeries);
len = len(2);
averagingNeighboursPerSide = 5; % on each side of an I-Frame, the averaging takes this many neighbours.
filteredTimeSeries = timeSeries; % This should create a deep copy of the input time series
%% Replace I-Frames with average of neighbours
for i = 1:len
    if (mod(i, iFramesInterval) == 1) %I-frames (if interval is 50, for example) are 1, 51, 101, ...
        n = 0;
        newValue = 0;
        for j = (i-1):-1:(i - averagingNeighboursPerSide)
            if j < 1
                break
            end
            newValue = newValue + filteredTimeSeries(j);
            n = n + 1;
        end
        for j = (i+1):(i + averagingNeighboursPerSide)
            if j > len
                break
            end
            newValue = newValue + filteredTimeSeries(j);
            n = n + 1;
        end
        filteredTimeSeries(i) = newValue/n;
    end
end
end

