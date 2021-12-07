% PUt here the IDs of mice you want to analyse
mouseID = { '177','184','446','453','457','235','239','255','339','342','349'};
% Put here the experimental conditions you want to analyse
experiment = {'baseline','ICA'};
% put here the length of each experimental condition (add 1 in the
% beginning for indexing)
explgt = [1 5 3];

%%
if exist('CAIM','var');clear CAIM;end
if exist('PCA','var');clear PCA;end
for i = 1:length(mouseID)    
    for j = 1:length(experiment)  
        %%
        % here the pathname is constructed were the data is
        pathname = ['/media/2Photon/Nicola/baselineday1/' mouseID{i} '/' experiment{j} '/*/'];
        % Here the analysis script is called. Go in there to choose your
        % analysis
        [CAIMtemp,PCAtemp] = CombBeCa(pathname,mouseID{i});
        % Here the index for the pooling variable (CAIM) is defined
        if j == 1
            ind = sum(explgt(1:j)):sum(explgt(2:j+1));
            ind = ind(length(ind)-length(CAIMtemp)+1:length(ind));
        elseif j ==2
            ind = sum(explgt(1:j)):sum(explgt(1:j))+length(CAIMtemp)-1;
        end

        CAIM(ind,i) = CAIMtemp;
        PCA(ind,i) = PCAtemp;
    end
end

 save(['/media/2Photon/Nicola/AnalysisFINAL/BigFatSummary.mat'],'CAIM');
%save(['/media/2Photon/Nicola/Analisi2020/BigFatPCA.mat'],'PCA');
 


