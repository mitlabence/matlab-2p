%% Include necessary packages
CAIMAN_PATH = '/home/bemitlas/Code/Matlab/external/CaImAn/';
CELLSORT_PATH = '/home/bemitlas/Code/Matlab/external/CellSort/';
NORMCORRE_PATH = '/home/bemitlas/Code/Matlab/external/Normcorre/';
EXTERN_PATH =  '/home/bemitlas/Code/Matlab/external/other/';
FILEHANDLING_PATH = '/home/bemitlas/Code/Matlab/FileHandling/';

addpath(genpath(CAIMAN_PATH)); %CaImAn replacing ca_source_extraction
addpath(genpath(CELLSORT_PATH));
addpath(genpath(NORMCORRE_PATH));
addpath(genpath(EXTERN_PATH));
addpath(genpath(FILEHANDLING_PATH));

%% Set files and parameters
DATA_PATH = '/media/bemitlas/raw-data/AG-Wenzel/Group/Test/T386_MatlabTest/'; 
FILE_ND2 = [DATA_PATH 'T386_20211202_green.nd2'];
FILE_ABF = [DATA_PATH '21d02000.abf'];
FILE_LABVIEW = [DATA_PATH 'T386.021221.1105.txt'];
FILE_TIME = [DATA_PATH 'T386.021221.1105time.txt'];


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
