


load(['/media/2Photon/Nicola/AnalysisFINAL/BigFatSummary.mat'])
experiment = {'Baseline'};
%%

% The coloums of CAIM refer to individual mice that were pooled. Set the following value to
% decide which mice should be included:
numice = [1:4];
% The rows refer to experimental conditions the were pooled. Decide here
% which experiments you want to pool
numexp = [5];
%%
ID = cell(0,length(mouseID),5);
if exist('CAIM','var');clear CAIM;end
if exist('PCA','var');clear PCA;end

for j = 1:length(experiment)
    Behave.mouseID = [];
    Behave.numrounds = [];
    Behave.runtime = [];
    Behave.pupilsize = [];   
    
    mouseIDgc = 1;
    mouseIDbulk = 1;
    
    %%
    for i = 1:length(mouseID)
        %%
        if ispc
         pathname = ['/media/2Photon/Nicola/rasterplot/' mouseID{i} '/' experiment{j} '/*/'];
        else
        pathname = ['/media/2Photon/Nicola/rasterplot/' mouseID{i} '/' experiment{j} '/*/'];
        end
        
          [files,samedate,prefiles] = FindDataSets(pathname);
        %[CAIMtemp,PCAtemp,caim,scn,behave] = CombBeCa(pathname,mouseID{i});

        %%
        mouseIDtemp = cell(size(CAIMtemp,1),4);
        mouseIDtemp(:,1) = mouseID(i);
        mouseIDtemp(:,2) = num2cell(behave.expID);
        mouseIDtemp(:,3) = num2cell(max(behave.expID));
        mouseIDtemp(:,4) = num2cell(i);
        mouseIDtemp(:,5) = experiment(j);
        
        Behave.mouseID = [Behave.mouseID; mouseIDtemp];
        Behave.numrounds = [Behave.numrounds; behave.numrounds];
        Behave.runtime = [Behave.runtime; behave.runtime];
        Behave.pupilsize = cat(3,Behave.pupilsize,behave.pupilsize);
        
        
        
        kk = size(mouseIDtemp,1);
        if i == 1
            k = size(ID,1)+1;
            ID = [ID; cell(kk,length(mouseID),5)];
%         elseif kk > size(ID,2)-k
%             ID = [ID; cell(kk-(size(ID,2)-k),length(mouseID),5)];
        end
        
        ID(k+cell2mat(mouseIDtemp(:,2))-1,i,:) = mouseIDtemp;
        CAIM(k:k+mouseIDtemp{end,3}-1,i) = CAIMtemp;
        PCA(k:k+mouseIDtemp{end,3}-1,i) = PCAtemp;
    end

end
%%
% if ispc
%     save(['C:\Users\Admin\Desktop\BigFatCluster.mat'],'CAIM','ID')
% else
%     save(['/media/2Photon/Martin Pofahl/BigFatCluster.mat'],'CAIM','ID')
%     save(['/media/2Photon/Martin Pofahl/BigFatPCA.mat'],'PCA');%,'-v7.3')
% end

