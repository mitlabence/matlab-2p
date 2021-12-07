clear;

%% load file

%TODO: check hard-coded values, and path to video and coordinates xlsx!

CAIMAN_PATH = 'D:\Software\CaImAn\CaImAn-MATLAB\'; %give with \ at the end

% test on small data
%DATA_PATH = 'D:\PhD\Matlab Scripts\testvid_caiman\'; %give with \ at the end
%file_name = 'movie_moco.tif';
%roi_file_name = 'moco_xlsx.xlsx';
%file_path = strcat(DATA_PATH, file_name);
%roi_file_path = strcat(DATA_PATH, roi_file_name);

% test on large file
DATA_PATH = 'D:\PhD\Matlab Scripts\testvid_large_caiman\';
file_name = 'T301_tmev_d1.270820.1209.tif';
% roi_file_name = 'Results.xlsx'; %TODO: create an xlsx file, or implement
% csv reading function...

addpath(genpath(strcat(CAIMAN_PATH, 'utilities'))); %add utilities library from CaImAn
%addpath(genpath('utilities')); %fixme: which version of utilities is this? Why does the newest one in CAIMAN folder not work?



centroids= []; %put in ROI centroid variable name
%centroids= coordinates; %TODO: find out what format this should be...
%already have csv file!
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=1500;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(file_path, sframe, num2read);

%Y = movie; %replace bigread2() if movie is already defined ("after
%bigread2" part of file)
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% Set parameters

K = 30;                                           % number of components to be found
tau = 4;                                          % std of gaussian kernel (size of neuron) 
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 1;                                  % merging threshold

% CNMFSetParams(): the fresh version (from CaImAn)
options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','ellipse','dist',3,...      % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );

% options.ROIList=xlsread(ROIListName);
% options.ROI_List=xlsread(roi_file_path); %FIXME: ROI_list should be in a parameter struct P passed in initialize_components!
%options.ROIList=[];

% FIXME: ROI_list takes (Y, X) formatted points, fiji exports (X, Y)...
% roi_opts.ROI_list=xlsread(roi_file_path);
roi_opts.ROI_list = readCoordinateList(roi_file_path, 'inverted');

%% Data pre-processing

[P,Y] = preprocess_data(Y,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options, roi_opts);  % initialize

% display centers of found components
Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

%% manually refine components (optional)
refine_components = false;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
    
%% update spatial components
Yr = reshape(Y,d,T);
clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain, bin],P,options);

%% update temporal components
[C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

%%
% display_merging = 1; % flag for displaying merging example
% if display_merging
%     i = 1; randi(length(merged_ROIs));
%     ln = length(merged_ROIs{i});
%     figure;
%         set(gcf,'Position',[300,300,(ln+2)*300,300]);
%         for j = 1:ln
%             subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2)); 
%                 title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
%         end
%         subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
%                 title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
%         subplot(1,ln+2,ln+2);
%             plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
%             hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
%             title('Temporal Components','fontsize',16,'fontweight','bold')
%         drawnow;
% end

%% repeat
[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,[Am, b],P,options);
[C2,f2,P,S2] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);

% inferred: spatial weighting on the ROI pixels, unmixing, background substraction, and denoising
% filtered: spatial weighting on the ROI pixels, unmixing, background substraction
% raw: uniform spatial weighting on the ROI pixels (with threshold to remove the very low energy pixels), background substraction
% df: df signal
% F: baseline
% dfof: dfof
 extractControl.baselineRatio=0.25;%the higher the more conservative, meaning the threshold becomes elevated%
[ inferred, filtered, raw ] = signalExtraction(Yr,A2,C2,b2,f2,d1,d2,extractControl);

inferSpikeControl.method='foopsi_derivative';
inferSpikeControl.stdThr=2.5;%the higher the more conservative, meaning the threshold becomes elevated%
inferSpikeControl.stdThr2=2.5;%the higher the more conservative, meaning the threshold becomes elevated%
inferSpikeControl.lowpassCutoff=1;            % practically a smoothing (the smaller the value the smoother; low pass in Hz (set to -1 if not performing low pass filter)
inferSpikeControl.frameRate=10;           % frame rate in Hz
inferSpikeControl.dynamicThr=1; %if this is zero, the whole thing only depends on one parameter, namely stdThr  
inferSpikeControl.SNR=max(C2,[],2)./cell2mat(P.neuron_sn);
inferSpikeControl.SNR0=12; %the higher the more conservative, meaning the threshold becomes elevated%
inferSpikeControl.multiplier=0.25;%the higher the more conservative%
inferSpikeControl.multiplier2=0.25;%the higher the more conservative%

spike = inferSpike( inferred.dfof, S2, inferSpikeControl);
figure; imagesc(spike)
colormap('gray');

plotControl2.frameRate=10;
plotControl2.normalization=1;   plotControl2.sep=1.2;   plotControl2.displayLabel=1;  plotControl2.rollingView=0;
plotControl2.plotInferred=1;    plotControl2.plotFiltered=1;   plotControl2.plotRaw=0;  plotControl2.plotSpike=1;

plotActivityTraceSpike(inferred.dfof(1:5,1:1000), filtered.dfof(1:5,1:1000), raw.dfof(1:5,1:1000), spike(1:5,1:1000), plotControl2);

% %% do some plotting
% 
% [A_or,C_or,S_or,P] = order_ROIs(A2,C2,S2,P); % order components
% K_m = size(C_or,1);
% [C_df,~,S_df] = extract_DF_F(Yr,[A_or,b2],[C_or;f2],S_or,K_m+1); % extract DF/F values (optional)
% 
% contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
% figure;
% [Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),contour_threshold,1); % contour plot of spatial footprints
% pause; 
% %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
% %% display components
% 
% plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options)
% 
% %% make movie
% 
% make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)
