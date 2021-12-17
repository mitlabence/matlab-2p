%% Include necessary packages

CAIMAN_PATH = 'D:\Software\CaImAn\CaImAn-MATLAB\';
CELLSORT_PATH = 'D:\Software\CellSort\CellSort\';
NORMCORRE_PATH = 'D:\Software\NoRMCorre\NoRMCorre\';
OTHER_EXTERNAL_PATH = 'D:\PhD\Matlab Ultron\codes_matlab\external\other\';
MARTIN_FUNCTIONS_PATH = 'D:\PhD\Matlab Ultron\codes_matlab\external\Martin Code\'; 
FILEHANDLING_PATH = 'D:\PhD\Matlab Ultron\codes_matlab\FileHandling\';
PROCESSING_PATH = 'D:\PhD\Matlab Ultron\codes_matlab\Processing\';

%IMPORTANT: See CNMFSetParams.m and @CNMF\CNMF.m for fields of CNMFParams
%and CNMF object; see NoRMCorreSetParms and
%@MotionCorrection\MotionCorrection.m for NormCorre fields!


%Long-term to-do list:
%TODO: memory map might be useful sometimes: https://de.mathworks.com/help/matlab/import_export/overview-of-memory-mapping.html
%TODO: it would be less chaotic to use CNMF and MotionCorrection objects
%(composite? mixin? subclass?) and custom classes (belt, LFP). Merging
%motioncorrection and CaImAn is complicated, not worth the clarity it 
%brings... in CaiMaN Python, they are already merged!
%TODO: add export-import capability to be able to switch between Python and
%Matlab seamlessly. Staying within CNMF framework is important also for
%this reason!

addpath(genpath(CAIMAN_PATH)); %CaImAn replacing ca_source_extraction
addpath(genpath(CELLSORT_PATH));
addpath(genpath(NORMCORRE_PATH));
addpath(genpath(OTHER_EXTERNAL_PATH));
addpath(genpath(MARTIN_FUNCTIONS_PATH));
addpath(genpath(FILEHANDLING_PATH));
addpath(genpath(PROCESSING_PATH));

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

K = 400;                                        % number of components to be found
tau = 8;                                        % std of gaussian kernel (radius of neuron in pixel)
p = 2;                                          % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                % merging threshold
should_refine_manually = false;                                 % Manually refine components

options_caiman = CNMFSetParms(...
    'split_data',0,...
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,...
    'spatial_method','regularized'...           %'constrained',...
    );



%options for non-caiman and non-normcorre functions.
%TODO: to avoid missing parameters, define a function like CNMFSetParms,
%and call it to supply default parameters!
options_other.filename = nd2_file_path;
options_other.sframe = 1;
options_other.numchan = 1; %FIXME: not used yet. Only using single channel nd2 readout!
options_other.crop_moco = [0 0 0 0]; %used to be [20 0 0 0]
options_other.crop_caiman = [30 30 30 30];
options_other.num2read = [];
%options.numchan  = numchan;
%options.sframe   = sframe;
%options.num2read = num2read;
%options.crop = crop_caiman;
bright_spikes = [];

%define options for normcorre
options_moco = NoRMCorreSetParms('d1',size(nd2_data,1),'d2',size(nd2_data,2),...
    'grid_size',[32,32],'mot_uf',4,'bin_width',200,...
    'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);

%% Read nd2 file
[nd2_data, options_other] = nd2ReadWithOptions(options_other);

%% Perform ripple noise reduction and movement correction

nd2_data = rippleRemovalWithOptions(nd2_data, bright_spikes, options_other);
[nd2_data,shifts_g,template,options_moco,col_shift] = normcorre_batch(nd2_data,options_moco);

%Save to tiff after motion correction and ripple noise removal (optional)
%writeToTiff(nd2_data, DATA_PATH, OUTPUT_NAME);


%% Create CaImAn object (instance of CMNF class)
%CNM.fit (see demo_script_class.m in CaImAn) does all the steps in one!

%This allows cleaner pipeline
CNM = CNMF;
CNM.optionsSet(options_caiman);
CNM.filename = options_other.filename; %just for documentation
CNM.fr = 15;    %15 Hz framerate. TODO: precise framerate is stored in nd2? Use nd2finfo? in preprocessBeltCaimInplace, frequency is printed.
%after motion correction, drop a band on the edge of the video, i.e. crop
CNM.loadArray(cropArray(nd2_data, options_other.crop_caiman));

%% Processing

CNM.preprocess;             % preprocessing (compute some quantities)
CNM.initComponents(K);      % initialization
CNM.plotCenters()           % plot center of ROIs detected during initialization

%TODO: optionally manually refine here!

CNM.updateSpatial();        % update spatial components
CNM.updateTemporal(0);      % update temporal components (do not deconvolve at this point)

%TODO: optionally component classification here. ML toolbox needed! 
%CNM.evaluateComponents();   % evaluate spatial components based on their correlation with the data
%CNM.CNNClassifier('')       % evaluate spatial components with the CNN classifier
%CNM.eventExceptionality();  % evaluate traces
%CNM.keepComponents();       % keep the components that are above certain thresholds

CNM.merge();
CNM.displayMerging();

%repeat processing
CNM.updateSpatial();
CNM.updateTemporal(); % This time use p in options of CNMF object
CNM.extractDFF();           % DETREND(!) and extract DF/F values.
%TODO: do we need deconvolve() at some point?
%TODO: detrend vs extract?

%do some plotting
figure;
CNM.plotContours();
CNM.plotComponentsGUI();     % display all components
%TODO: in csem.m, which does roughly the same steps as this section, there
%is a mysterious Y_r that replaces Y. It originates in plot_components_GUI.
%what is it?

%Apparently, backgrounds (CNM.b) is a flattened 2d matrix, no temporal
%component. (in the case of 30-30-30-30 crop factor, 452*452 = 204304 x 1
%array.
%One entry that is weird is Df, which is empty at this point!

save([DATA_PATH OUTPUT_FILE_NAME 'Ca.mat'], 'CNM')
save([DATA_PATH OUTPUT_FILE_NAME 'moco.mat'], 'options_moco'); % save the options

%% Custom analysis

%TODO: end of csem.m: sort_components, divcells and canorm perform 
%non-trivial things, but are not documented. Decrypt them!
[cID, thresh] = sort_components(CNM); %seems not to change CMN object
%caim = divcells(CNM, cID); %FIXME: this would change CMN! Return a raw
%data instead of adding "raw" field to CMNF object.
%[S_norm, S_bin] = canorm(CNM);

%TODO: Check preprocessBeltCaimInplace: returns CNMF object?
%TODO: check BeltToSCN: returns CNMF objects?
%FIXME: belt has a modified tsscn (577 points, matched to nikon frames), 
%but everything else is kept untouched(?)...

%[belt, CNM, nikon_time_stamps, labview_time_stamps] = openImagingSession(DATA_PATH,OUTPUT_FILE_NAME, CNM);

%[belt, caim] = processBeltCaimInplace(belt, CNM, nikon_time_stamps, labview_time_stamps); 

%[belt,scn] = BeltToSCN(caim,belt); %scn is the belt in scn time frame
%scn_reduced = rmfield(scn, 'pupilsize'); %all other fields should have matching dimensions, hence writable to a single csv
%writetable(struct2table(scn_reduced), 'D:\PhD\Data\T386_MatlabTest\scn_reduced.csv');
