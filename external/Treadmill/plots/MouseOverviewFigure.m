mouseID = { 'prove' };
experiment = {'342a'};
%%
for i = 1:length(mouseID)    
    for j = 1:length(experiment)  
        %%
        pathname = ['/media/2Photon/Nicola/AnalysisFINAL/' mouseID{i} '/' experiment{j} '/*/'];
               
        [files,samedate,prefiles] = FindDataSets(pathname);
        
        %% here is where the read out happens
        for k = 1:size(samedate,1)
            %%
            disp(['load ' prefiles(k).name])
            load(prefiles(k).name);           
            
            SummaryPlot([files(k).folder '/'],files(k).name(1:end-6),caim,belt,scn)
            if ~isfield(caim,'FOV')               
                tiffile = dir(['/media/2Photon/Nicola/**/*' files(k).name(1:end-9) '*.tif']);
                if length(tiffile) > 1
                    disp('Multiple tiff files present. Delete obsolete tiff files')
                else
                    options.crop = [0 0 0 0]; % [70 0 10 5] for M227  % [70 0 0 5] for M234
                    options.sframe = 1;
                    options.num2read = 200;
                    [Data,Datared] = readdata([tiffile(1).folder '/' tiffile(1).name],options);
                    caim.FOV = mat2gray(mean(Data,3));
                    caim.FOV(:,:,2) = mat2gray(mean(Datared,3));
                    save([files(samedate{i,1}).folder '/preprocessed/' files(samedate{i,1}(1)).name(1:end-6) 'pro.mat'],'caim','belt','scn')
                end
            end
        end
    end
end