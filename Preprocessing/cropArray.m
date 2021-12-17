function cropped_array = cropArray(input_array, crop_xy)
%CROPARRAY given a 2d or 3d array (x-y-t, for example), crop in the x-y
%plane, throwing away: along x, the first crop_xy(1) and the last
%crop_xy(2), along y: the first crop_xy(3) and the last crop_xy(4) pixels.
%   Input:
%       input_array: 3d array with 3d dimension that will remain uncropped
%       crop_xy: [x1, x2, y1, y2] array, see functionality above.
%   Output:
%       cropped_array: the cropped array
cropped_array = input_array(crop_xy(1)+1:end-crop_xy(2), crop_xy(3)+1:end-crop_xy(4), :);
end

