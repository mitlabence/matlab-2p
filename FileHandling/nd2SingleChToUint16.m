function im_ch1 = nd2SingleChToUint16(filename,sframe,num2read)

FIRST_IMAGE_FRAME_IN_FILE_SYSTEM = 8; %see nd2info(filename).file_structure for location of imaging frames

%%
tic
finfo = nd2finfo(filename);

if nargin<2 || isempty(sframe)
    sframe = 1;
end
if nargin<3 || isempty(num2read)
    num2read = finfo.img_seq_count;
end
%%
im_ch1 = zeros(finfo.img_width, finfo.img_height, num2read,'uint16');

fid = fopen(filename, 'r');
filepos = find(strncmp('ImageDataSeq|0!',{finfo.file_structure(:).nameAttribute}, 14),1);
filepos = filepos+sframe-1;
fseek(fid,finfo.file_structure(filepos).dataStartPos, 'bof');
%%
tic
% Image extracted from ND2 has image width defined by its first dimension.
if finfo.padding_style == 1
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(2, :);
        fseek(fid, 2, 'cof');
    end
else
    if finfo.ch_count == 2
        %%
        disp('2 channels detected in nd2 file. Channel 2 lost!');
        for kk = 1:num2read
            if strcmp(finfo.file_structure(filepos).nameAttribute(1:10),'CustomData')
                filepos = filepos + 1;
            end
            fseek(fid,finfo.file_structure(filepos).dataStartPos+8, 'bof');
            for ii = 1: finfo.img_height
                temp = fread(fid, finfo.ch_count * finfo.img_width, '*uint16');
                temp = reshape(temp,[finfo.ch_count finfo.img_width]);
                im_ch1(:, ii,kk) = temp(1, :);
            end
            filepos = filepos + 1;
        end
    elseif finfo.ch_count == 1
        %%
        frameDataLength = finfo.file_structure(FIRST_IMAGE_FRAME_IN_FILE_SYSTEM).dataLength; %normal length of data corresponding to an imaging frame
        disp('1 channel found. Started reading out nd2 file.');
        frame_number = 1;
        while frame_number < (num2read+1) %equivalent to original for frame_number 1:num2read
            %TODO: see finfo.file_structure. The data we want are
            %those corresponding to 'ImageDataSeq[]!', which start at
            %position 8. In between, there are 'CustomDataSeq[] !'
            currentDataLength = finfo.file_structure(filepos).dataLength;
            if(currentDataLength ~= frameDataLength)
                filepos = filepos + 1;
                %Problem is, we should only increase frame number if we
                %actually write into the im_ch1 data! If this executes,
                %it means there was a 'CustomDataSeq' frame or whatnot, not
                %a real frame.
                continue
            end
            fseek(fid,finfo.file_structure(filepos).dataStartPos+8, 'bof');
            for ii = 1: finfo.img_height
                temp = fread(fid, finfo.ch_count * finfo.img_width, '*uint16');
                temp = reshape(temp,[finfo.ch_count finfo.img_width]);
                im_ch1(:, ii,frame_number) = temp(1, :);               
            end
            filepos = filepos + 1;
            frame_number = frame_number + 1;
        end
    end
end

fclose(fid);

im_ch1 = permute(im_ch1, [2 1 3]);

disp(['reading complete image data used ', sprintf('%0.2f', toc), ' seconds'])
end