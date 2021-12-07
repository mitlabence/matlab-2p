function outputTimeSeries = MagicalManipulation(derivativeTimeSeries)
%MAGICALMANIPULATION Summary of this function goes here
%   Detailed explanation goes here
outputTimeSeries = derivativeTimeSeries;
for i=1:(length(outputTimeSeries)-1)
    outputTimeSeries(i) = outputTimeSeries(i) - outputTimeSeries(i+1);
end
end

