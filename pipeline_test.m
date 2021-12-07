%% Include necessary packages
CAIMAN_PATH = '/home/bemitlas/Code/Matlab/external/CaImAn/';
CELLSORT_PATH = '/home/bemitlas/Code/Matlab/external/CellSort/';
NORMCORRE_PATH = '/home/bemitlas/Code/Matlab/external/Normcorre/';
OTHER_EXTERNAL_PATH = '/home/bemitlas/Code/Matlab/external/other/';

%FIXME: writing to Tiff - some frames are bad. When does it go wrong?
%During ripple removal, or during writing to tiff?

%FIXME: nd2read causes buggy frames like 257! Mostly 252-253 frames apart

addpath(genpath(CAIMAN_PATH)); %CaImAn replacing ca_source_extraction
addpath(genpath(CELLSORT_PATH));
addpath(genpath(NORMCORRE_PATH));
addpath(genpath(OTHER_EXTERNAL_PATH));

%% Set files and parameters
DATA_PATH = '/media/2Photon/AG Wenzel/tmev/T329/T329_tmev_d4/';
FILE_ND2 = [DATA_PATH 'T329_tmev_d4_05122020_002.nd2'];
FILE_ABF = [DATA_PATH '20d05004.abf'];
FILE_LABVIEW = [DATA_PATH 'T329_tmev_d4.051220.1224.txt'];
FILE_TIME = [DATA_PATH 'T329_tmev_d4.051220.1224time.txt'];

OUTPUT_NAME = 'test_out'; %will be a tif file


sframe = 1;
numchan = 1;
crop = [20 0 0 0];
num2read = [];

amplitudeThreshold = [10.8 12.8];
win = [40 40];

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

%% Read nd2
%For now, only 1 channel. 
%TODO: readdata() is there to handle different file types! For now, .nd2 is
%enough. But later .mat, .tif should be included too!

Data = nd2SingleChToUint16(FILE_ND2,sframe,num2read);


%% Ripple noise, motion correction
[Data,~] = RippleNoiseRemoval(Data,amplitudeThreshold(1),win(1),0,[]);

%motion correction parameters
options_nonrigid = NoRMCorreSetParms('d1',size(Data,1),'d2',size(Data,2),...
        'grid_size',[32,32],'mot_uf',4,'bin_width',200,...
        'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
[Data,~,~,~] = normcorre_batch(Data,options_nonrigid);

%% Save to .tif
writeToTiff(Data, DATA_PATH, OUTPUT_NAME);