
% mouseID = {'M103', 'M155','M158','M194','M195','M224','M226','M227','M229','M234'};
mouseID = { '184','177'};
% experiment = {'Baseline','Cues','Odor','LateBase','Airpuff'};
experiment = {'Baseline'};

i = 1;j = 1;

if ispc
    pathname = ['Z:/media/2Photon/Nicola/plot/' mouseID{i} '/' experiment{j} '/*/'];
    addpath(genpath('C:\Users\Admin\Documents\MATLAB\ca_source_extraction'))
else
    pathname = ['Z:/media/2Photon/Nicola/plot/'  mouseID{i} '/' experiment{j} '/*/'];
    addpath(genpath('/home/martin/Documents/MATLAB/ca_source_extraction'))
end

%%

 if exist('CAIM','var');clear CAIM;end
 if exist('PCA','var');clear PCA;end
 
 for j = 1:length(experiment)  
    for j = 1:length(experiment)     Behave.mouseID = [];
    Behave.numrounds = [];
    Behave.runtime = [];
    Behave.pupilsize = [];
    
    Fire.meanfire = [];
    Fire.frfld = [];
    Fire.fireprob = [];
    Fire.amplitude = [];
    Fire.meanamplitude = [];
    Fire.mouseID = [];
   Fire.placecode = [];
    Fire.numcells = [];
    Fire.cclust = [];
%    Fire.shufflep = [];

%    Network.netprob = [];
%    Network.netresp = [];
%    Network.netsum = [];
%    Network.netnum = [];
%    Network.netfreq = [];
%    Network.netoccur = [];
%    Network.netfirecorr =[];
%    Network.mouseID = [];
%    Network.replay = [];

%    Bulk.bulkbase = [];
%    Bulk.bulkstd = [];
%    Bulk.bulkratio = [];
%    Bulk.bulkspeed = [];
%    Bulk.bulknet = [];
%    Bulk.netxbulk = [];
 %   Bulk.bulkind = [];
 %   Bulk.netfld = [];
 %   Bulk.mouseID = [];
 %   Bulk.inout = [];
    
%    Response.runonset = RespClust;   
 %   Response.spcresp =  RespClust;
 %   Response.airpuff =  RespClust;
 %   Response.network =  RespClust;
 %   Response.soundl =   RespClust;
 %   Response.soundr =   RespClust;
 %   Response.stripes =  RespClust;
    
  %  mouseIDgc = 1;
  %  mouseIDbulk = 1;
    
    %%
    for i = 1:length(mouseID)
        %%
        if ispc
              pathname = ['/media/2Photon/Nicola/Analisi2020/' mouseID{i} '/' experiment{j} '/*/'];
        else
            pathname = ['/media/2Photon/imaging Negar & Martin/' mouseID{i} '/Analysis/' experiment{j} '/'];
        end
        
        [behave,fire,network,bulk,caim,scn,response,PCAtemp] = CombBeCa(pathname,mouseID{i});
        
%         fire.frfld = fire.frfld(sum(fire.frfld,2)~=0,:);
        %%
        mouseIDtemp = cell(size(fire.frfld,1),4);
        mouseIDtemp(:,1) = mouseID(i);
        mouseIDtemp(:,2) = num2cell(behave.expID);
        mouseIDtemp(:,3) = num2cell(max(behave.expID));
        mouseIDtemp(:,4) = num2cell(i);
        mouseIDtemp(:,5) = experiment(j);
        mouseIDtempgc = mouseIDtemp; 
        mouseIDtempbulk = mouseIDtemp;
        
        Behave.mouseID = [Behave.mouseID; mouseIDtemp];
        Behave.numrounds = [Behave.numrounds; behave.numrounds];
        Behave.runtime = [Behave.runtime; behave.runtime];
        Behave.pupilsize = cat(3,Behave.pupilsize,behave.pupilsize);
        
        if isfield(caim,'A')                 
                     
            mouseIDtempgc(:,4) = num2cell(mouseIDgc); 
            mouseIDgc = mouseIDgc + 1;
            Fire.mouseID = [Fire.mouseID; mouseIDtempgc];
            Network.mouseID = [Network.mouseID; mouseIDtempgc];
            
            Fire.numcells =  [Fire.numcells;fire.numcells];
            Fire.meanfire =  [Fire.meanfire; fire.meanfire];
            Fire.frfld =     [Fire.frfld; fire.frfld];
            Fire.fireprob =  [Fire.fireprob; fire.fireprob];
            Fire.amplitude = [Fire.amplitude; fire.amplitude];
            Fire.meanamplitude = [Fire.meanamplitude; fire.meanamplitude];           
            Fire.placecode = [Fire.placecode; fire.placecode];
            Fire.cclust =    [Fire.cclust; fire.cclust];
            Fire.shufflep = [Fire.shufflep; fire.shufflep];
            
            Network.netprob = [Network.netprob network.netprob];
            Network.netnum = [Network.netnum network.netnum];
            Network.netresp = cat(3,Network.netresp,network.netresp);
            Network.netsum = [Network.netsum network.netsum];
            Network.netfreq = [Network.netfreq; network.netfreq];
            Network.netoccur = [Network.netoccur; network.netoccur];
            Network.netfirecorr = [Network.netfirecorr network.netfirecorr];
            Network.replay = [Network.replay; network.replay];
                               
        end
        
        if ~isempty(bulk.baseline)                      
            
            mouseIDtempbulk(:,4) = num2cell(mouseIDbulk);
            mouseIDbulk = mouseIDbulk + 1;
            Bulk.mouseID = [Bulk.mouseID; mouseIDtempbulk];
            
            Bulk.bulkbase = [Bulk.bulkbase; bulk.baseline];
            Bulk.bulkstd = [Bulk.bulkstd; bulk.std];
            Bulk.bulkratio = [Bulk.bulkratio; bulk.ratio];
            Bulk.bulkspeed = [Bulk.bulkspeed; bulk.speedcorr];
            if isfield(caim,'A')
                Bulk.bulknet = [Bulk.bulknet; bulk.netcorr];
                Bulk.netxbulk = [Bulk.netxbulk; bulk.netxbulk];
                Bulk.bulkind = [Bulk.bulkind; bulk.indcorr];
                Bulk.netfld = [Bulk.netfld; bulk.netfld];
                if isfield(response,'airpuff')
                    Bulk.inout = [Bulk.inout; [response.network.inout(1:size(response.airpuff.inout),:) response.airpuff.inout]];
                end
            end    
        end
        
        if isfield(response,'airpuff')
            Response.airpuff = RespClust(Response.airpuff,response.airpuff,[isfield(caim,'A') ~isempty(bulk.baseline)],mouseIDtempgc,mouseIDtempbulk);
        end
        
        if isfield(response,'runonset')
            Response.runonset = RespClust(Response.runonset,response.runonset,[isfield(caim,'A') ~isempty(bulk.baseline)],mouseIDtempgc,mouseIDtempbulk);
        end
        
        if isfield(response,'spcresp') && ~isempty(response.spcresp.times)
            Response.spcresp= RespClust(Response.spcresp,response.spcresp,[isfield(caim,'A') ~isempty(bulk.baseline)],mouseIDtempgc,mouseIDtempbulk);
        end
        
        if isfield(response,'network')
            Response.network = RespClust(Response.network,response.network,[isfield(caim,'A') ~isempty(bulk.baseline)],mouseIDtempgc,mouseIDtempbulk);
        end
        
        if isfield(response,'soundl')
            Response.soundl= RespClust(Response.soundl,response.soundl,[isfield(caim,'A') ~isempty(bulk.baseline)],mouseIDtempgc,mouseIDtempbulk);
        end
        
        if isfield(response,'soundr')
            Response.soundr= RespClust(Response.soundr,response.soundr,[isfield(caim,'A') ~isempty(bulk.baseline)],mouseIDtempgc,mouseIDtempbulk);
        end
        
        if isfield(response,'stripes')
            Response.stripes= RespClust(Response.stripes,response.stripes,[isfield(caim,'A') ~isempty(bulk.baseline)],mouseIDtempgc,mouseIDtempbulk);
        end
        
        
        kk = size(mouseIDtemp,1);
        if i == 1
            k = size(ID,1)+1;
            ID = [ID; cell(kk,length(mouseID),5)];
%         elseif kk > size(ID,2)-k
%             ID = [ID; cell(kk-(size(ID,2)-k),length(mouseID),5)];
        end
        
        ID(k+cell2mat(mouseIDtemp(:,2))-1,i,:) = mouseIDtemp;
        CAIM(k:k+mouseIDtemp{end,3}-1,i) = caim;
        PCA(k:k+mouseIDtemp{end,3}-1,i) = PCAtemp;
    end
    
    %%    
    if ispc
        save(['C:\Users\Admin\Desktop\' experiment{j} '.mat'],'Behave','Fire','Network','Bulk','Response')
    else
        save(['/media/2Photon/Martin Pofahl/' experiment{j} '.mat'],'Behave','Fire','Network','Bulk','Response')
    end

end
%%
if ispc
    save(['C:\Users\Admin\Desktop\BigFatCluster.mat'],'CAIM','ID')
else
    save(['/media/2Photon/Martin Pofahl/BigFatCluster.mat'],'CAIM','ID')
    save(['/media/2Photon/Martin Pofahl/BigFatPCA.mat'],'PCA');%,'-v7.3')
end

toc