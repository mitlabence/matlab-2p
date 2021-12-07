addpath(genpath('/home/fabian/Documents/MATLAB/ca_source_extraction'));
addpath(genpath('/media/2Photon/Fabian/MATLAB/NoRMCorre-master'));
addpath(genpath('/media/2Photon/Fabian/MATLAB/sorter'));
addpath(genpath('/media/2Photon/Fabian/Current/'));


PathName = '/media/2Photon/Nicola/177-184/184/184.151118/second10/';
files = dir([PathName '.tif']);
tic;
%% Readout and cropping
num2read = [];
sframe = 1;
numchan = 2;
crop = [20 0 0 0]; %for ripple noise
% crop = [0 0 0 0]; %for NNMF

%% Ripple noise removal

amplitudeThreshold = [10.8 12.8];
win = [40 40];
% RippleNoiseRemoval(Data,amplitudeThreshold(1),win(1),1,[]);
% RippleNoiseRemoval(Datared,amplitudeThreshold(2),win(2),1,[]);


%% Set parameters

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
options.crop = crop;

%%
for i = 1:length(files)
    close all
    disp(['Working on ' files(i).name])
    nam = [files(i).folder '/' files(i).name];
    [Data,Datared] = readdata(nam,options);
    %% Ripple noise
    [Data,~] = RippleNoiseRemoval(Data,amplitudeThreshold(1),win(1),0,[]);
    [Datared,~] = RippleNoiseRemoval(Datared,amplitudeThreshold(2),win(2),0,[]);
    
    %% Movement Correction
    options_nonrigid = NoRMCorreSetParms('d1',size(Data,1),'d2',size(Data,2),...
        'grid_size',[32,32],'mot_uf',4,'bin_width',200,...
        'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
    
    [Datared,shifts,template,options_nonrigid] = normcorre_batch(Datared,options_nonrigid); %channel used to perform the movement correction

    % Apply shifts to the red channel (channel to which movement correction withh be applied)
    
    Data = apply_shifts(Data,shifts,options_nonrigid);

 
    %% NNMF
    
    crop = [30 30 30 30]; % NNMF
    caim = csem(Data(crop(1)+1:end-crop(2),crop(3)+1:end-crop(4),:),K,p,refine,options);
    if size(Data,3)<200;int = size(Datam,3);else; int = 200;end
    caim.FOV = mean(Data(:,:,1:int),3);
    caim.FOV(:,:,2) = mean(Datared(:,:,1:int),3);
    
    
    %% Bulk signal
    % Threshold should be high enough that binary plot (on the left)
    % resembles average plot in the middle
    %          bulk.cropCA1 = [0 350 0 400]; % CA1 input
    %     bulk.threshCA1 = 0.2;
    %     % %     [bulk.traceCA1,bulk.templateCA1] = bulkdata(Datared,Data,bulk.threshCA1,bulk.cropCA1,2,0);
    %     bulk.cropMEC = [0 0 0 0]; % MEC input
    %     bulk.threshMEC = [];
    %     [bulk.traceMEC,bulk.templateMEC] = bulkdata(Datared,Data,bulk.threshMEC,bulk.cropMEC,2,0);
    % %     bulk.cropMC = [200 80 0 0]; % MC layer input
    %     bulk.threshMC = 0.2;
    %     [bulk.traceMC,bulk.templateMC] = bulkdata(Datared,Data,bulk.threshMC,bulk.cropMC,2,0);
    %     caim.bulk = bulk;
    
    %% Save
    
    TiffData(Data,Datared,files(i).folder,files(i).name);
    
    a = find(nam == '.');
    save([nam(1:a(end)-1) 'Ca.mat'],'caim')
    clear Data Datared
end
toc