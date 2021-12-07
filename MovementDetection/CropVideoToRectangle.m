function CropVideoToRectangle(input_fpath, cropRectangle, output_fpath)
%CROPVIDEOTORECTANGLE Given an input file input_fpath and a rectangle,
%saves a cropped video to output_fpath.
%   input_fpath: full path and filename of video to crop (video will not be
%   modified!)
%   rectangle: a list of [x, y, w, h]
%   output_fpath: full path including output file name.

x = cropRectangle(1);
y = cropRectangle(2);
width = cropRectangle(3);
height = cropRectangle(4);

system(sprintf('ffmpeg -i %s -filter:v "crop=%i:%i:%i:%i" %s', input_fpath, width, height, x, y, output_fpath));

end

