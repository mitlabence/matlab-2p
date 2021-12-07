function caim = csem(Y,K,p,refine,options)%% load file 
%This function seems to be getting single cells using CaImAn. Necessary to
%perform this function before passing the result to BeltToSCN? There,
%matching might be done using this caim properties, that are not present
%before (after Ripple Noise removal and motion correction).
%% Data pre-processing
[d1,d2,T] = size(Y);
d = d1*d2;                   % total number of pixels
options.d1 = d1;
options.d2 = d2;
tau = options.gSig;
[P,Y] = preprocess_data(Y,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

% display centers of found components
Cn =  reshape(P.sn,d1,d2); %correlation_image(Y); %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
figure;imagesc(Cn);
    axis equal; axis tight; hold all;
    scatter(center(:,2),center(:,1),'mo');
    title('Center of ROIs found from initialization algorithm');
    drawnow;

%% manually refine components (optional)
refine_components = false;%refine;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end
    
%% update spatial components
Yr = reshape(Y,d,T);
clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

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
[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,[Am,b],P,options);
[C2,f2,P,S2] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);

%% do some plotting

[A_or,C_or,S_or,P] = order_ROIs(A2,C2,S2,P); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,[A_or,b2],[C_or;f2],K_m+1); % extract DF/F values (optional)

%contour_threshold = 0.95;                       % amount of energy used for each component to construct contour plot
% figure;
% [Coor,json_file] = plot_contours(A_or,reshape(P.sn,d1,d2),options,0); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

%% display components
[Y_r,Df] = plot_components_GUI_modifiedByNico(Yr,A_or,C_or,b2,f2,Cn,options);


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