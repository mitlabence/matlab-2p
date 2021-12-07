function [im_ch1, im_ch2, im_ch3] = nd2read(filename,sframe,num2read)
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
im_ch2 = zeros(finfo.img_width, finfo.img_height, num2read,'uint16');
im_ch3 = zeros(finfo.img_width, finfo.img_height, num2read,'uint16');
if finfo.ch_count == 4
  im_ch4 = zeros(finfo.img_width, finfo.img_height, 'uint16');
end

fid = fopen(filename, 'r');
filepos = find(strncmp('ImageDataSeq|0!',{finfo.file_structure(:).nameAttribute}, 14),1);
filepos = filepos+sframe-1;
fseek(fid,finfo.file_structure(filepos).dataStartPos, 'bof');
%%
tic
% Image extracted from ND2 has image width defined by its first dimension.
if finfo.padding_style == 1
  if finfo.ch_count == 4
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
%         im_ch3(:, ii) = temp(1, :);
        im_ch1(:, ii) = temp(2, :);
        im_ch2(:, ii) = temp(3, :);
%         im_ch4(:, ii) = temp(4, :);
        fseek(fid, 2, 'cof');
    end
  else
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch3(:, ii) = temp(1, :);
        im_ch1(:, ii) = temp(2, :);
        im_ch2(:, ii) = temp(3, :);
        fseek(fid, 2, 'cof');
    end
  end
else
    if finfo.ch_count == 4
        for ii = 1: finfo.img_height
            temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
              [finfo.ch_count finfo.img_width]);
            im_ch1(:, ii) = temp(1, :);
            im_ch2(:, ii) = temp(2, :);
            im_ch3(:, ii) = temp(3, :);
            im_ch4(:, ii) = temp(4, :);
        end
    elseif finfo.ch_count == 3
        for ii = 1: finfo.img_height
            temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
              [finfo.ch_count finfo.img_width]);
            im_ch1(:, ii) = temp(1, :);
            im_ch2(:, ii) = temp(2, :);
            im_ch3(:, ii) = temp(3, :);
        end 
    elseif finfo.ch_count == 2
        %%
        for kk = 1:num2read
            if strcmp(finfo.file_structure(filepos).nameAttribute(1:10),'CustomData')
                filepos = filepos + 1;
            end
            fseek(fid,finfo.file_structure(filepos).dataStartPos+8, 'bof');
            for ii = 1: finfo.img_height
                temp = fread(fid, finfo.ch_count * finfo.img_width, '*uint16');
                temp = reshape(temp,[finfo.ch_count finfo.img_width]);
                im_ch1(:, ii,kk) = temp(1, :);
                im_ch2(:, ii,kk) = temp(2, :);
            end
            filepos = filepos + 1;
        end
    elseif finfo.ch_count == 1
        %%
        for kk = 1:num2read
            fseek(fid,finfo.file_structure(filepos).dataStartPos+8, 'bof');
            for ii = 1: finfo.img_height
                temp = fread(fid, finfo.ch_count * finfo.img_width, '*uint16');
                temp = reshape(temp,[finfo.ch_count finfo.img_width]);
                im_ch1(:, ii,kk) = temp(1, :);               
            end
            filepos = filepos + 1;
        end
    end
end

fclose(fid);

im_ch1 = permute(im_ch1, [2 1 3]);
im_ch2 = permute(im_ch2, [2 1 3]);
im_ch3 = permute(im_ch3, [2 1 3]);
% if finfo.ch_count == 4
%     im_ch4 = permute(im_ch4, [2 1]);
% end
% if any(strcmpi(varargin, 'use_ch4'))
%   im_ch3 = im_ch4;
% end
  

disp(['reading complete image data used ', sprintf('%0.2f', toc), ' seconds'])
end