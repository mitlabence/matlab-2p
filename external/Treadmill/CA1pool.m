mouseID = {'M177','M184'};%
experiment = {'Baseline'};
clear CAIM PCA
k = 1;
for i = 1:length(mouseID)
    
    for j = 1:length(experiment)  
        %%
        if ispc
            pathname = ['Z:\Martin Pofahl\CA1\' mouseID{i} '\Analysis\' experiment{j} '\'];
        else
            pathname = ['/media/2Photon/Martin Pofahl/CA1/' mouseID{i} '/Analysis/' experiment{j} '/'];
        end
        
        [behave,fire,~,~,caim,~,~,PCAtemp] = CombBeCa(pathname,mouseID{i});
        mouseIDtemp = cell(size(fire.frfld,1),4);
        mouseIDtemp(:,1) = mouseID(i);
        mouseIDtemp(:,2) = num2cell(behave.expID);
        mouseIDtemp(:,3) = num2cell(max(behave.expID));
        mouseIDtemp(:,4) = num2cell(i);
        mouseIDtemp(:,5) = experiment(j); 
        
        ID(k+cell2mat(mouseIDtemp(:,2))-1,i,:) = mouseIDtemp;
        CAIM(k:k+mouseIDtemp{end,3}-1,i) = caim;
        PCA(k:k+mouseIDtemp{end,3}-1,i) = PCAtemp;
        %% here is where the read out happens
%         [files,samedate,prefiles] = FindDataSets(pathname);      
%         for k = 1:size(samedate,1)        
%             disp(['load ' prefiles(k).name])
%             load(prefiles(k).name);           
%             
%             SummaryPlot(pathname,files(samedate{k,1}(1)).name(1:end-6),caim,belt,scn)
%             PCAprint(caim,scn,pathname,files(samedate{k,1}(1)).name(1:end-6),'pca');
%         end
    end
end

if ispc
    save(['Z:\Martin Pofahl\CA1BigFatCluster.mat'],'CAIM')
else
    CA1 = CAIM; PCA1 = PCA;
    save(['/media/2Photon/Martin Pofahl/CA1/BigFatCluster.mat'],'CA1','PCA1')
end