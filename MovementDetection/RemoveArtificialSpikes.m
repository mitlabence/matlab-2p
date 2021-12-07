function l_filtered = RemoveArtificialSpikes(vectors_list)
%REMOVEARTIFICIALSPIKES Removes the spikes appearing at every 50th frame
%(the I-frames), which is a consequence of how compressed videos are
%stored. NOTE: originally intended for use on time series from
%estimateOpticalFlow of cropped video frames, it does not work, because the
%optical flow algorithm introduces weird delays. For example, in one video,
%the first 52 spikes coincide with the original I-frames, then there is a 1
%frame delay between I-frames and the shown spikes in the optical flow time
%series. This goes on for 106 I-frames, then the delay increases to 2 for
%the next 105, again for 106, then again for 105, then again for 106 ...
%In another video, the first 51 I-frames show spikes, then 1 delay in the
%optical flow for the next 106, then one more for the next 105 ...
%Extremely weird and stressful! Decided to implement a magical filter
%instead.
%   vectors_list: the per-frame magnitude of motion vectors of a movie

DROP_PERCENT = 0.02; % drop 2% of highest values (as every 50th frame is aberrant, and they are higher than signal)


s = size(vectors_list);
N = s(1); %number of frames

sorted_l = sort(vectors_list); % increasing order
cut_off = sorted_l(floor((1.0-DROP_PERCENT)*N)); % roughly every 50th frame, i.e. 2% shows the artefacts. Get corresponding cut-off value.

clear s;
clear sorted_l;

l_filtered = vectors_list;

for i = 1:N
    if l_filtered(i) >= cut_off
        % below, select j = previous element (cyclical, i.e. first
        % element's previous neighbour is last element), k = next element
        if i == 1
            j = N;
            k = 2;
        else
            if i == N
                j = i-1;
                k = 1;
            else
                j = i-1;
                k = i+1;
            end
        end
        l_filtered(i) = ( l_filtered(j) + l_filtered(k) ) / 2.0; %create average of neighbors
    end
end

end
