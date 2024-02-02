function [im_ch1, im_ch2, im_ch3] = nd2ToUint16(filename, sframe, channels_to_read, num2read)
% Supports reading out up to 3 channels (channel 4 seems to be a special
% case in nd2read.m)
% channels_to_read: 1-indexing list, for example [1, 2] to read both
% channels, [1] to read only channel 1, [2] only channel 2...

FIRST_IMAGE_FRAME_IN_FILE_SYSTEM = 8; %see nd2info(filename).file_structure for location of imaging frames

%%
tic
finfo = nd2finfo(filename);
n_channels = finfo.ch_count;

if nargin<2 || isempty(sframe)
    sframe = 1;
end

if nargin<3 || isempty(channels_to_read)
    % read out all channels
    channels_to_read = 1:n_channels;
else
    % user specified which channels to read out
    % take only existing channels from the list the user specified
    channels_to_read = sort(channels_to_read(channels_to_read <= n_channels));
end

if nargin<4 || isempty(num2read)
    num2read = finfo.img_seq_count;
end

%%
if ismember(1, channels_to_read)
    im_ch1 = zeros(finfo.img_width, finfo.img_height, num2read,'uint16');
else
    im_ch1 = [];
end
if ismember(2, channels_to_read)
    im_ch2 = zeros(finfo.img_width, finfo.img_height, num2read,'uint16');
else
    im_ch2 = [];
end
if ismember(3, channels_to_read)
    im_ch3 = zeros(finfo.img_width, finfo.img_height, num2read,'uint16');
else
    im_ch3 = [];
end

if isempty(channels_to_read)  % empty list of channels explicitely specified
    return;
end

fid = fopen(filename, 'r');
filepos = find(strncmp('ImageDataSeq|0!',{finfo.file_structure(:).nameAttribute}, 14),1);
filepos = filepos+sframe-1;
fseek(fid,finfo.file_structure(filepos).dataStartPos, 'bof');
%%
tic
% Image extracted from ND2 has image width defined by its first dimension.
if finfo.padding_style == 1
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, n_channels * finfo.img_width, '*uint16'),...
          [n_channels finfo.img_width]);
        if ismember(1, channels_to_read)
            im_ch3(:, ii) = temp(1, :);
        end
        if ismember(2, channels_to_read)
            im_ch1(:, ii) = temp(2, :);
        end
        if ismember(3, channels_to_read)
            im_ch2(:, ii) = temp(3, :);
        end
        fseek(fid, 2, 'cof');
    end
else
    if n_channels == 3
        for ii = 1: finfo.img_height
            temp = reshape(fread(fid, n_channels * finfo.img_width, '*uint16'),...
              [n_channels finfo.img_width]);
            if ismember(1, channels_to_read)
                im_ch1(:, ii) = temp(1, :);
            end
            if ismember(2, channels_to_read)
                im_ch2(:, ii) = temp(2, :);
            end
            if ismember(3, channels_to_read)
                im_ch3(:, ii) = temp(3, :);
            end
        end 
    elseif n_channels == 2
        %%
        for kk = 1:num2read
            if strcmp(finfo.file_structure(filepos).nameAttribute(1:10),'CustomData')
                filepos = filepos + 1;
            end
            fseek(fid,finfo.file_structure(filepos).dataStartPos+8, 'bof');
            for ii = 1: finfo.img_height
                temp = fread(fid, n_channels * finfo.img_width, '*uint16');
                temp = reshape(temp,[n_channels finfo.img_width]);
                if ismember(1, channels_to_read)
                    im_ch1(:, ii,kk) = temp(1, :);
                end
                if ismember(2, channels_to_read)
                    im_ch2(:, ii,kk) = temp(2, :);
                end
            end
            filepos = filepos + 1;
        end
    elseif n_channels == 1
        %%
        for kk = 1:num2read
            fseek(fid,finfo.file_structure(filepos).dataStartPos+8, 'bof');
            for ii = 1: finfo.img_height
                temp = fread(fid, n_channels * finfo.img_width, '*uint16');
                temp = reshape(temp,[n_channels finfo.img_width]);
                im_ch1(:, ii,kk) = temp(1, :);               
            end
            filepos = filepos + 1;
        end
    end
end

fclose(fid);

if ismember(1, channels_to_read)
    im_ch1 = permute(im_ch1, [2 1 3]);
end
if ismember(2, channels_to_read)
    im_ch2 = permute(im_ch2, [2 1 3]);
end
if ismember(3, channels_to_read)
    im_ch3 = permute(im_ch3, [2 1 3]);
end
disp(['reading complete image data used ', sprintf('%0.2f', toc), ' seconds'])
end