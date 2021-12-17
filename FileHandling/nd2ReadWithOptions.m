function [im_ch1, options] = nd2ReadWithOptions(options)
%ND2READWITHOPTIONS Based on nd2SingleChToUint16, taking options
%structure for better code uniformity (avoid using different parameters for
%consequent steps, by supplying one option structure and possibly modifying
%it when options are missing). At the end, the user can look at the options
%and retrace what was done with which options.
%   Input:
%       options.filename: complete file name and path of nd2 file. Example:
%           'C:\Nikon\example_recording.nd2'
%       options.sframe: starting frame to read out. Default: 1
%       options.num2read: number of frames to read out. Default: length of
%           file.
%   Output:
%       im_ch1: d1 x d2 x n_frames Uint16 array.
%       options: the (possibly modified) input options structure.

if ~isfield(options,'filename') || isempty(options.filename) 
    error('Error in nd2ReadWithOptions: no filename supported in options!'); 
end
if ~isfield(options,'sframe') || isempty(options.sframe)
    disp('nd2ReadWithOptions: Supplying default parameter sframe: 1');
    options.sframe = 1; %starting frame is 1 (i.e. from beginning) by default
end
if ~isfield(options,'num2read') || isempty(options.num2read)
    disp('nd2ReadWithOptions: Supplying default parameter num2read: []');
    options.num2read = [];
end

im_ch1 = nd2SingleChToUint16(options.filename, options.sframe, options.num2read);
end

