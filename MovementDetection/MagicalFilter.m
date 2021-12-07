function filteredTimeSeries = MagicalFilter(timeSeries)
%MAGICALFILTER This filter uses some magic to remove the I-frame artifacts
%from a time series.
%   Detailed explanation goes here
%% Create return time series and first derivative of time series
IframePeriod = 50; % the period between I-frames. In HIKVision movies, this is
% 50, i.e. frames 1, 51, 101, ... are I-frames.
filteredTimeSeries = timeSeries;
derivatives = DerivativeOfTimeSeries(timeSeries); % the spikes at i turn into positive spikes at i and negative at i+1.
magicalTimeSeries = MagicalManipulation(derivatives); % experimentally promising manipulation of derivatives to find peaks
% basically, y(i) = y(i) - y(i+1). At the outlier position i, the y(i+1) is
% a large negative value while y(i) is a large positive value. Even if the
% outlier is not the largest value and a simple maximum-search would not 
%find it, subtracting the next value (the large negative value) will
%increase it, while mostly for "natural" high values, the time series looks
%more continuous, hence there won't likely be a sudden huge change in first
%derivative.
%Even more basically: calculate the negative second derivative. It is
%positive for maxima, and we expect the outliers to be a more outstanding
%peak than the natural ones, even if the outlier peak is surrounded by
%peaks that are larger.

averagingNeighboursPerSide = 2; % use this many neighbours on each side to replace outlier with average

%% Try to find outlier in each window
for i=1:IframePeriod:length(timeSeries)
   upperLimit = min((i+IframePeriod-1), length(timeSeries));
   window =  magicalTimeSeries(i:upperLimit); % choose the interval such that only 1
   % I-frame is in it (ideally).
   %find maximum and corresponding index
   [~, max_index] = max(window);
   %if y(i) is the outlier, then y(i)-y(i+1) will be large, hopefully the
   %maximum
   outlierIndex = i+max_index-1;
   %create average of neighbours
    n = 0;
    newValue = 0;
    for j = (outlierIndex-1):-1:(outlierIndex - averagingNeighboursPerSide)
        if j < 1
            break
        end
        newValue = newValue + filteredTimeSeries(j);
        n = n + 1;
    end
    for j = (outlierIndex + 1):(outlierIndex + averagingNeighboursPerSide)
        if j > length(filteredTimeSeries)
            break
        end
        newValue = newValue + filteredTimeSeries(j);
        n = n + 1;
    end
    filteredTimeSeries(outlierIndex) = newValue/n;
end

%% Remove I-frame artifacts from second compression (every 250th frame), which is apparently undistorted by optical flow...
for i=1:250:length(filteredTimeSeries)
    n = 0;
    newValue = 0;
    for j = (i-1):-1:(i - averagingNeighboursPerSide)
        if j < 1
            break
        end
        newValue = newValue + filteredTimeSeries(j);
        n = n + 1;
    end
    for j = (i + 1):(i + averagingNeighboursPerSide)
        if j > length(filteredTimeSeries)
            break
        end
        newValue = newValue + filteredTimeSeries(j);
        n = n + 1;
    end
    filteredTimeSeries(i) = newValue/n;
end
end

