function derivatives = DerivativeOfTimeSeries(timeSeries)
%DERIVATIVEOFTIMESERIES Summary of this function goes here
%   Calculate the derivative as follows:
%   dy/dx(i) = (y(i) - y(i-1))/dx. dx currently is taken as 1 for the sake
%   of simplicity.
%   The very first element will use y(0) := 0, i.e. dy/dx(1) = y(1).
derivatives = zeros(size(timeSeries));
derivatives(1) = timeSeries(1);

for i=2:length(derivatives)
   derivatives(i) = timeSeries(i) - timeSeries(i-1); 
end
end

