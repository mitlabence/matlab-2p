function [filtered_data,bright_spikes] = rippleRemovalWithOptions(nd2_data, bright_spikes, options)
%RIPPLEREMOVALWITHOPTIONS Performs ripple noise removal with supplied
%options as a structure. Uses the external function RippleNoiseRemoval.m.
%   Input:
%       nd2_data: 3d uint16 array (e.g. from nd2ReadWithOptions.m)
%       bright_spikes: list of bright spikes. Supply empty array to let
%           algorithm run from beginning.
%       options: options structure, of which the following are used:
%           options.rnr_ampl_thr: amplitude threshold
%           options.rnr_win: window size for algorithm
%           options.rnr_plt: boolean, whether to plot
%   Oputput:
%       filtered_data: the resulting rippple noise filtered data
%       bright_spikes: list of bright spikes, see algorithm.

if ~isfield(options,'rnr_ampl_thr') || isempty(options.filename) 
    disp('Supplying default parameter rnr_ampl_thr: 10.8');
    options.rnr_ampl_thr = 10.8;
end
if ~isfield(options,'rnr_win') || isempty(options.filename) 
    disp('Supplying default parameter rnr_win: 40');
    options.rnr_win = 40;
end
if ~isfield(options,'rnr_plt') || isempty(options.filename) 
    disp('Parameter rnr_plt has not been provided; will not plot! Setting rnr_plt to false.');
    options.rnr_plt = false;
end

[filtered_data, bright_spikes] = RippleNoiseRemoval(nd2_data,options.rnr_ampl_thr,options.rnr_win,options.rnr_plt,bright_spikes);
%As next step, use motion correction (normcorre_batch.m, for example).
end

