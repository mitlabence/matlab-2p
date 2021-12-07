function coordinates = readCoordinateList(fpath,mode)
%READCOORDINATELIST Summary of this function goes here
%   Read the supplied file (fpath) of two columns, X and Y, and return the
%   list either as-is or inverted (1st column Y, 2nd column X).
%   fpath: string path of xlsx file.
%   mode: 'normal' or 'inverted'. Default is 'normal', meaning first column
%           returned as X, second as Y.
if isempty(mode)
   mode = 'normal'; 
end
coords_list = xlsread(fpath);
list_shape = size(coords_list);
if list_shape(2) ~= 2
   error('Error in readCoordinateList: xlsx contains more (or less) than 2 columns (contains %d)!', list_shape(2));
end
if strcmp(mode, 'inverted')
   x_ind = 2;
   y_ind = 1;
else
   x_ind = 1;
   y_ind = 2;
end

first_col = coords_list(:, x_ind);
second_col = coords_list(:, y_ind);

coordinates = [first_col, second_col];