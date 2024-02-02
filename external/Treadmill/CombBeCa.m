function [CAIM,PCA] = CombBeCa(pathname,mouseID)

[files,samedate,prefiles] = FindDataSets(pathname);
if isempty(files)
   CAIM = [];
   PCA = [];
    return
end
expID = cell2mat(samedate(:,2));
%%
for i = 1:size(samedate,1) 
    close all
    %% Load Data if already processed before
    if ~isempty(prefiles) && length(prefiles)>=i
        disp(['load ' prefiles(i).name])
        load(prefiles(i).name,'belt','caim','scn');       
    else
        %% Readout ca-Imaging & behaviour data and real time correction    
        [belt,caim] = ConnecCaim([files(samedate{i,1}).folder '/'],files,samedate{i,1}); 
        if isfield(caim, "caim")  % fix(?) bug where caim becomes a struct containing caim itself only
            caim = caim.caim;
        end
        [belt,scn] = BeltToSCN(caim,belt);   
    end
    
    %% Place and Speed coding analysis
    scn = placecor(caim,scn);
    %% Firing properties analysis
         caim.fireprop = firefreq(caim,scn);
    %% Network events analysis     
         caim.network = netevents(caim,scn); 
    %% Shuffle Analysis
 %    caim = bulkanalysis(caim,scn);
        scn = stimcor(belt,caim,scn);
  %   numit = 1000;
        runshuf = 1;
  %   caim.shuffle = ShuffleAna(belt,caim,scn,numit,runshuf);
  %   caim.netthresh = caim.shuffle.netthresh;
  %   caim.network = netevents(caim,scn); 
    %% PCA analysis
       runshuf = false;
      
        caim.PCAout = PCAana(caim,scn,'pca');
        caim.GLMout = GLMana(caim,scn,'pca',runshuf);
        caim.ICAout = PCAana(caim,scn,'ica');
 %    caim.GPFAout = PCAsectioned(caim,scn,'gpfa');
    %% Input Bulk Signal analysis
%      caim = bulkanalysis(caim,scn);
    %% Correlation bulk & network signal
%      caim =  NetCorr(caim,scn);    
    %% Correlation to stimuli
%      scn = stimcor(belt,caim,scn);
    %% Cluster analysis of network events
       caim.network = ClustAna(caim,scn);      
    %% Cell cluster analysis   
    
      scn.cclust = cellcluster(caim,scn,str2double(mouseID(2:end)));       

    %% Pooling struct output
    
    if isfield(scn,'cellID') && size(scn.spaccode,2)>1
        plcfield = NaN(size(scn.spaccode,1),size(scn.spaccode,2));
        scn.plcfield;plcfield(scn.cellID(:,1),:) = scn.plcfield(:,:,1);
    elseif isfield(caim,'Y')
        plcfield = NaN(size(caim.Y,1),150);
    else
        plcfield = [];
    end
    
    CAIM(expID(i)).behave.numrounds = max(scn.rounds);
    CAIM(expID(i)).behave.totdist = scn.totdist(end-1);
    CAIM(expID(i)).behave.meanspeed = 100*mean(scn.speed(scn.running==1));
    CAIM(expID(i)).behave.runtime = scn.tsscn(end)*sum(scn.running==1)/length(scn.tsscn)/1000;
    CAIM(expID(i)).behave.pupilsize = scn.pupilsize;
    CAIM(expID(i)).behave.running = logical(scn.running);
    CAIM(expID(i)).behave.speed = scn.speed;
    CAIM(expID(i)).behave.tsscn = scn.tsscn;
    
    if isfield(caim,'Y')
        CAIM(expID(i)).Y = [];
        CAIM(expID(i)).A = caim.A;
        CAIM(expID(i)).S = logical(caim.S_bin);
        CAIM(expID(i)).Cn = caim.Cn;
        CAIM(expID(i)).cclust = scn.cclust;
        CAIM(expID(i)).plcfield = plcfield;
        CAIM(expID(i)).network = caim.network;
        CAIM(expID(i)).speedcorr = scn.speedcorr;
        CAIM(expID(i)).fireprop = caim.fireprop;
        PCA(expID(i)).PCAout = caim.PCAout;
        PCA(expID(i)).ICAout = caim.ICAout;
        PCA(expID(i)).GPFAout = [];
        PCA(expID(i)).GLMout = caim.GLMout;
    else
        
        CAIM(expID(i)).Y = [];
        CAIM(expID(i)).A = [];
        CAIM(expID(i)).S = [];
        CAIM(expID(i)).Cn = [];
        CAIM(expID(i)).cclust = [];
        CAIM(expID(i)).plcfield = [];
        CAIM(expID(i)).network = [];
        CAIM(expID(i)).speedcorr = [];
        CAIM(expID(i)).fireprop = [];
        PCA(expID(i)).PCAout = [];
        PCA(expID(i)).ICAout = [];
        PCA(expID(i)).GPFAout = [];
        PCA(expID(i)).GLMout = [];
    end
    
    if isfield(caim,'bulk')
        CAIM(expID(i)).bulk = caim.bulk;
    else
        CAIM(expID(i)).bulk = [];
    end   
    %% save preproceesed data
     if ispc
         save([files(samedate{i,1}).folder '\preprocessed\' files(samedate{i,1}(1)).name(1:end-6) 'pro.mat'],'caim','belt','scn')
     else
         save([files(samedate{i,1}).folder '/preprocessed/' files(samedate{i,1}(1)).name(1:end-6) 'pro.mat'],'caim','belt','scn')
    end
    
    %% Summary Plot and .pdf
  %    SummaryPlot(pathname,files(samedate{i,1}(1)).name(1:end-6),caim,belt,scn);
%      PCAprint(caim,scn,pathname,[files(samedate{i,1}(1)).name(1:end-6) '_PCA'],'pca');
end
end
 
    
    

