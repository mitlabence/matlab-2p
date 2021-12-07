%% Include necessary packages
CAIMAN_PATH = 'D:\Software\CaImAn\CaImAn-MATLAB\';
CELLSORT_PATH = 'D:\Software\CellSort\CellSort\';
NORMCORRE_PATH = 'D:\Software\NoRMCorre\NoRMCorre\';
OTHER_EXTERNAL_PATH = 'D:\PhD\Matlab Ultron\codes_matlab\external\other\';
MARTIN_FUNCTIONS_PATH = 'D:\PhD\Matlab Ultron\codes_matlab\external\Martin Code\'; 
FILEHANDLING_PATH = 'D:\PhD\Matlab Ultron\codes_matlab\FileHandling\';

%FIXME: writing to Tiff - some frames are bad. When does it go wrong?
%During ripple removal, or during writing to tiff?

%FIXME: nd2read causes buggy frames like 257! Mostly 252-253 frames apart

%TODO: memory map might be useful sometimes: https://de.mathworks.com/help/matlab/import_export/overview-of-memory-mapping.html

addpath(genpath(CAIMAN_PATH)); %CaImAn replacing ca_source_extraction
addpath(genpath(CELLSORT_PATH));
addpath(genpath(NORMCORRE_PATH));
addpath(genpath(OTHER_EXTERNAL_PATH));
addpath(genpath(MARTIN_FUNCTIONS_PATH));
addpath(genpath(FILEHANDLING_PATH));

%% Set files and parameters
DATA_PATH = 'D:\PhD\Data\T386_MatlabTest\'; %folder in which the data is located; should end with \
ND2_FNAME = 'T386_20211202_green.nd2';

FILE_ND2 = [DATA_PATH 'T386_20211202_green.nd2']; %Complete path of nd2
FILE_NIKDATA = [DATA_PATH 'T386.021221.1105nik.txt']; %Complete path of experiment data exported from .nd2 file
FILE_ABF = [DATA_PATH '21d02000.abf']; %Abf file (LFP) complete path
FILE_LABVIEW = [DATA_PATH 'T386.021221.1105.txt']; %LabView-generated txt file with columns for velocity, distance etc.
FILE_TIME = [DATA_PATH 'T386.021221.1105time.txt']; %LabView-generated txt file with time stamps for matching with Nikon

% File name (without extension); if "xy", then "xyCa.mat" is output of Ca
% data after ripple noise removal, motion correction and single-cell
% labeling.
OUTPUT_FILE_NAME = 'T386.021221.1105'; 


nd2_file_path = [DATA_PATH ND2_FNAME];
sframe = 1;
numchan = 1;
crop_moco = [0 0 0 0]; %used to be [20 0 0 0]
crop_caiman = [30 30 30 30];
num2read = [];

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

options.numchan  = numchan;
options.sframe   = sframe;
options.num2read = num2read;
options.crop = crop_caiman;


%% Read nd2 file
nd2_data = nd2SingleChToUint16(nd2_file_path,sframe,num2read);

%% Perform ripple noise reduction and movement correction
nd2_data = correctNoiseMovement(nd2_data, sframe, numchan, crop_moco, num2read);

%writeToTiff(nd2_data, DATA_PATH, OUTPUT_NAME);

%% Single cell stuff
caim = csem(nd2_data(crop_caiman(1)+1:end-crop_caiman(2), crop_caiman(3)+1:end-crop_caiman(4),:), K, p, refine, options);
save([DATA_PATH OUTPUT_FILE_NAME 'Ca.mat'], 'caim')


[belt,caim] = readcaim(DATA_PATH, OUTPUT_FILE_NAME);%FIXME: belt has a modified tsscn
%(577 points, matched to nikon frames), but everything else is kept
%untouched(?)...
%FIXME: readcaim() is a terrible function! It makes lots of changes not
%obvious from the function name... readxy() should not ever change the
%read-out data!

[belt,scn] = BeltToSCN(caim,belt); %scn is the belt in scn time frame

scn_reduced = rmfield(scn, 'pupilsize'); %all other fields should have matching dimensions, hence writable to a single csv
writetable(struct2table(scn_reduced), 'D:\PhD\Data\T386_MatlabTest\scn_reduced.csv');