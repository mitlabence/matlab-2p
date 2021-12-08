function caim = csem(Y,K,p,refine,options)
% This function seems to be getting single cells using CaImAn. Necessary to
% perform this function before passing the result to BeltToSCN? There,
% matching might be done using this caim properties, that are not present
% before (after Ripple Noise removal and motion correction).
%   input:
%       Y: series of image frames (_PROBABLY_ uint16 entries) d1 x d2 x T
%       K: suggested number of neurons to extract (presumably not forced?)
%       p: order of autoregressive system in NNMF (y_t = a_0 + a_1 * y_t-1
%       + ... + a_p * y_t-p + A_t (A_t white noise)
%       refine: boolean; should manual refinement take place?
%   options    fine-tuning parameters (optional)
%       Used by:
%           * initialize_components (Caiman)
%           * manually_refine_components (Caiman)
%           * update_spatial_components (Caiman)
%           * update_temporal_components (Caiman)
%           * merge_components (Caiman)
%           * plot_components_GUI_modifiedByNico (custom, taken from Caiman)
%       options will be set as a field of caim.
%       
%       Fields: (same as in CNMFSetParms)
%         Defined in initialize_components:
%           options.init_method: method of initialization ('greedy','sparse_NMF','HALS')
%           options.gSig: half size of neurons to be found (default: [5,5])
%                           Also set as tau here, used by initialize_components
%                           for other parameters as default.
%           options.nIter: number of iterations for shape tuning (default 5)
%           options.gSiz: size of kernel (default 2*tau + 1)
%           options.ssub: spatial downsampling factor (default 1)
%           options.tsub: temporal downsampling factor (default 1)
%           options.nb: rank of background component (default 1)
%           options.save_memory: flag for processing data in chunks to save memory (default 0)
%           options.windowSiz: size of spatial window when computing the median (default 32 x 32)
%           options.chunkSiz: number of timesteps to be processed simultaneously if on save_memory mode (default: 100)
%           options.med_app: number of timesteps to be interleaved for fast (approximate) median calculation (default: 1, no approximation)
%           options.rem_prct: percentile to be removed before initialization (default: 20)
%           options.noise_norm: normalization by noise estimate prior to initialization (default: true)
%           options.noise_norm_prctile: minimum noise level (as percentile of P.sn) used in the normalization prior to initialization (default: 2)
%         Defined in manually_refine_components:
%           options.d1: first dimension of frame (set in this function)
%           options.d2: second dimension of frame (set in this function)
%           options.cont_threshold: ?? (documentation says it is used in
%               plot_contours.m, but it is used here, without explanation.
%         Defined in update_spatial_components.m:
%           options.interp
%           options.spatial_parallel
%           options.search_method
%         Defined in update_temporal_components.m: ...
%         Defined in merge_components.m: ...
%           See CaImAn CNMFSetParams.m
%   output:
%       caim: CNMF object, see @CNMF/CNMF.m
% Original author: Martin Pofahl(?) (copied from code collection of Nicola 
% Masala in 2021, original function name: csem)
% Refactoring: Bence Mitlasoczki, 2021

%% Data pre-processing
[d1,d2,T] = size(Y);
d = d1*d2; % total number of pixels per frame
options.d1 = d1;
options.d2 = d2;
tau = options.gSig;
[P,Y] = preprocess_data(Y,p);

%% Fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);

%% Display centers of found components
Cn =  reshape(P.sn,d1,d2); 

%correlation_image(Y); 
%max(Y,[],3); 
%std(Y,[],3); % image statistic (only for display purposes)

figure;
    imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
drawnow;

%% manually refine components (optional)
if refine
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
    
    
%% update spatial components
Y_flattened = reshape(Y,d,T); %flatten Y (d1 x d2 x T) to a 2D matrix (d x T).
clear Y; %FIXME: this does not remove Y from workspace!
[A,b,Cin] = update_spatial_components(Y_flattened,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S] = update_temporal_components(Y_flattened,A,b,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Y_flattened,A,b,C,f,P,S,options);

%%
display_merging = 0; % flag for displaying merging example
if and(display_merging, ~isempty(merged_ROIs))
    i = 1; %randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    figure;
        set(gcf,'Position',[300,300,(ln+2)*300,300]);
        for j = 1:ln
            subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2)); 
                title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
        end
        subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
                title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
        subplot(1,ln+2,ln+2);
            plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
            hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
            title('Temporal Components','fontsize',16,'fontweight','bold')
        drawnow;
end

%% repeat
P.p = p;    % restore AR value
[A2,b2,Cm] = update_spatial_components(Y_flattened,Cm,f,[Am,b],P,options);
[C2,f2,P,S2] = update_temporal_components(Y_flattened,A2,b2,Cm,f,P,options);

%% do some plotting

[A_or,C_or,S_or,P] = order_ROIs(A2,C2,S2,P); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Y_flattened,[A_or,b2],[C_or;f2],K_m+1); % extract DF/F values (optional)

%contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
% figure;
% [Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),options,0); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

%% display components
[Y_r,Df] = plot_components_GUI_modifiedByNico(Y_flattened,A_or,C_or,b2,f2,Cn,options);


%% save

caim.Y = Y_r;
caim.A = A_or;
caim.C = C_or;
caim.S = S_or;
caim.b = b2;
caim.f = f2;
caim.Cn = Cn;
caim.Df = Df;
caim.options = options;
[caim.cID,caim.thresh] = sort_components(caim);
caim  = divcells(caim,caim.cID);
[caim.S_norm,caim.S_bin] = canorm(caim);



end