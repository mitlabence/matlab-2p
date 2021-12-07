function [time_list, optical_flow_time_series] = GetOpticalFlow(video_path)
%GETOPTICALFLOW Given a (.avi, if possible) video, this function returns a
%1D time series where the points correspond to the average magnitude of
%motion vectors over the frame.
%   video_path: the complete path to the video to get the time series from.

vidReader = VideoReader(video_path);
opticFlow = opticalFlowHS;

N = vidReader.NumFrames;
T = vidReader.Duration;
dt = T/N; % dt is the time between two frames
optical_flow_time_series = zeros(N, 1);
time_list = (0:dt:(T-dt)); % create time points for data points

for i = 1:N
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);  
    flow = estimateFlow(opticFlow,frameGray); %FIXME: as the videos contain I-frames every 50th and 250th frame,
    % the delay on the flow-artifacts related to the I-frames should come
    % from this algorithm somehow!
    avg = mean(flow.Magnitude, 'all'); % average magnitude of flow vectors over frame
    optical_flow_time_series(i) = avg;
    if mod(i, 100) == 0
        disp(i);
    end
end


end

