function [filteredData,brightSpikes] = RippleNoiseRemoval(Data,amplitudeThreshold,win,plt,brightSpikes)
%% Filter out periodic ripple noise from a 2-Ph image
 

if plt
    figure('position',[36 56 1692 647])
    grayImage = Data(:,:,1);%-min(grayImage(:));

    % Compute the 2D fft to find ripple noise
    frequencyImage = fftshift(fft2(grayImage));

    if isempty(brightSpikes)
        % Take log magnitude 
        amplitudeImage = log(abs(frequencyImage));
        h = subplot(2,4,[1 5]); % 2 rows, 4 columns; 1 and 5 make up the left column
        imagesc(amplitudeImage)
        colormap(h,gray)
%         minValue = min(min(amplitudeImage));
%         maxValue = max(max(amplitudeImage));
%         imshow(amplitudeImage, [minValue maxValue]);
%         axis on;

        brightSpikes = amplitudeImage > amplitudeThreshold;

        %Create rectangle of window to aid manual refinement of selecting
        %components to filter out
        rectangle_filter_boundary = zeros(size(amplitudeImage));
        rectangle_filter_boundary(round(end/2-win):round(end/2+win), round(end/2-win)) = 1;
        rectangle_filter_boundary(round(end/2-win):round(end/2+win), round(end/2+win)) = 1;
        rectangle_filter_boundary(round(end/2-win), round(end/2-win):round(end/2+win)) = 1;
        rectangle_filter_boundary(round(end/2+win), round(end/2-win):round(end/2+win)) = 1;
        
        %FIXME: not sure what happens if brightSpikes[i,j] is not 0, and 1
        %is added! Hopefully the color does not change by much.
        subplot(2,4,2);imagesc(brightSpikes + rectangle_filter_boundary); %FFT spectrum (log abs value, thresholded) before filtering? Also plot rectangle of filtering
        %imshow(brightSpikes); 
        
        subplot(2,4,3);imagesc(amplitudeImage);
%         colormap(jet)
        % Remove the central DC spike, exclude everything from row 115 to 143
        brightSpikes(round(end/2-win):round(end/2+win),round(end/2-win):round(end/2+win)) = 0;
        subplot(2,4,6);
        imagesc(brightSpikes + rectangle_filter_boundary); % FFT spectrum (log abs value, thresholded) after filtering? Also plot rectangle of filtering
        subplot(2,4,7);
        hold on
        for i_row = 1:size(amplitudeImage, 1)
            plot(amplitudeImage(i_row,:));
        end
        yline(amplitudeThreshold);
        hold off
    end

    % Mask amplitude spectrum
    % frequencyImage = fftshift(fft2(grayImage));
    frequencyImage(brightSpikes) = 0;

    %Get filtered image using Inverse fast Fourier transform
    % filteredImage = ifft2(fftshift(frequencyImage));
    filteredImage = ifft2((frequencyImage));

    %here is the filtered image
    filteredData = abs(filteredImage);

    %visualize input image(grayImage) and output image(amplituImage3)
    %this part can be commented out for when the program is used in a loop

    subplot(2,4,4);
    imagesc(grayImage);
%     minValue = min(grayImage(:));
%     maxValue = max(grayImage(:));
%     imshow(grayImage, [minValue maxValue]);
    title('Input Image');
    subplot(2,4,8);
    imagesc(filteredData);
%     minValue = min(filteredData(:));
%     maxValue = max(filteredData(:));
%     imshow(filteredData, [minValue maxValue]);
    title('Filtered Image');
end

if ~plt && isempty(brightSpikes)
    filteredData = zeros(size(Data));
    parfor frame=1:size(Data,3)
        grayImage = Data(:,:,frame);%-min(grayImage(:));
        frequencyImage = fftshift(fft2(grayImage));
        amplitudeImage = log(abs(frequencyImage));
        brightSpikes = amplitudeImage > amplitudeThreshold; 
        brightSpikes(round(end/2-win):round(end/2+win),round(end/2-win):round(end/2+win)) = 0;
        frequencyImage(brightSpikes) = 0;
        filteredImage = ifft2((frequencyImage));
        filteredData(:,:,frame) = abs(filteredImage);
    end
    brightSpikes = [];
elseif ~plt && ~isempty(brightSpikes)   
    filteredData = zeros(size(Data));
    parfor frame=1:size(Data,3) 
        grayImage = Data(:,:,frame);%-min(grayImage(:));
        frequencyImage = fftshift(fft2(grayImage));
        frequencyImage(brightSpikes) = 0;
        filteredImage = ifft2((frequencyImage));
        filteredData(:,:,frame) = abs(filteredImage);
    end
    brightSpikes = [];
end



end
