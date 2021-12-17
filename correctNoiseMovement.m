function nd2_data = correctNoiseMovement(uint16_data, sframe, numchan, crop, num2read)
%CORRECTNOISEMOVEMENT This function combines the steps of ripple noise
%removal and motion correction.
%   Input:
%       uint16_data: 3d array of uint16 elements (stack of imaging frames,
%           e.g. from nd2SingleChToUint16.m
%       sframe: the skip_frame parameter ofa CNMF object. TODO: it seems to
%           not be used!
%       numchan: number of channels
% 
%FIXME: options is not used?!
%TODO: as part of the function unification process, change input to a CNMF
%object!

%% Ripple noise removal
amplitudeThreshold = [10.8 12.8];
win = [40 40];

[nd2_data, ~] = RippleNoiseRemoval(uint16_data,amplitudeThreshold(1),win(1),0,[]);

%% Motion correction

K = 400;                                        % number of components to be found
tau = 8;                                        % std of gaussian kernel (radius of neuron in pixel)
p = 2;                                          % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                % merging threshold
refine = false;                                 % Manually refine components


options = CNMFSetParms(...
    'split_data',0,...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,...
    'spatial_method','regularized'...           %'constrained',...
    );

%not used, for reading nd2
options.numchan  = numchan;
options.sframe   = sframe;
options.num2read = num2read;
options.crop = crop;


options_nonrigid = NoRMCorreSetParms('d1',size(nd2_data,1),'d2',size(nd2_data,2),...
        'grid_size',[32,32],'mot_uf',4,'bin_width',200,...
        'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
[nd2_data,~,~,~] = normcorre_batch(nd2_data,options_nonrigid); %[nd2_data,shifts,template,options_nonrigid, col_shift]
% if a second channel should be corrected too: Data = apply_shifts(Data,shifts,options_nonrigid);

%% nd2_data is returned after processing
end

