function options = setOptionsForPython(filename)
%SETOPTIONSFORPYTHON This function makes it easy in Python to prepare a
% struct options for nd2ReadWithOptions.
%   Detailed explanation goes here

%In python, folder separators are / (slash). Change them to Matlab
%separators
options.filename = strrep(filename, '/', '\');
end

