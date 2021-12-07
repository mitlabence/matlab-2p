function cclust = cellcluster(caim,scn,expIDin)

cclustID.expID =             1;     % 1 mouse ID
cclustID.meanf =             2;     % 2 mean frequency of cell firing
                                    % 3 mean frequency of cell firing during run
                                    % 4 mean frequency of cell firing during rest
cclustID.plcfld =            5;     % 5 placefield center
cclustID.plcfldp =           6;     % 6 placefield p value
cclustID.plclength =         7;     % 7 length of place field
cclustID.plcvctang =         8;     % 8 place coding vector angle
cclustID.plcvct =            9;     % 9 place coding vector length
cclustID.plcvctp =          10;     % 10 place coding vector p value
cclustID.plcfirst =         11;     % 11 first round with place field activity
cclustID.speedcorr =        12;     % 12 logical if cell is speedmodulated
cclustID.speedcorrp1 =      13;     % 13 p value of pearson correlation
cclustID.speedcorrp2 =      14;     % 14 p value of shuffle analysis
cclustID.spatialIN =        15;     % 15 Skaggs spatial information, baseline corrected
cclustID.spatialINp =       16;     % 16 Skaggs spatial information, baseline corrected


cclustID.netperc =          17;     % 17 percentage of network related firing
cclustID.netnon =           18;     % 18 percentage of NON network related firing
cclustID.netprob =          19;     % 19 netprob
cclustID.clust =            20;     % 20 cluster ID
cclustID.coactsum =         21;     % 21 number of coactive cells
cclustID.coactmean =        22;     % 22 mean number of coactive events 
cclustID.corrpossum =       23;     % 23 number of sig pos correlated cells
cclustID.corrposmean =      24;     % 24 mean of sig pos correlations
cclustID.corrnegsum =       25;     % 25 number of sig neg correlated cells
cclustID.corrnegmean =      26;     % 26 mean of sig neg correlations
cclustID.coactplacesum =    27;     % 27 number of coactive place cells
cclustID.coactplacemean =   28;     % 28 mean of coactive place cell - events
cclustID.corrplacepossum =  29;     % 29 number of sig pos correlated place cells
cclustID.corrplaceposmean = 30;     % 30 mean of sig pos place correlations
cclustID.corrplacenegsum =  31;     % 31 number of sig neg correlated place cells
cclustID.corrplacenegmean = 32;     % 32 mean of sig neg place correlations

cclustID.nevent =           33;     % 33 # stim-response events
cclustID.nstim =            34;     % 34 # responded stimuli
cclustID.stimang =          35;     % 35 place preference angle of stim response
cclustID.stimleng =         36;     % 36 length of stim-ppv
cclustID.coactstimsum =     37;     % 37 number of coactive stim cells
cclustID.coactstimmean =    38;     % 38 mean of coactive place cell - events
cclustID.corrstimpossum =   39;     % 39 number of sig pos correlated place cells
cclustID.corrstimnegsum =   40;     % 40 number of sig neg correlated place cells
cclustID.corrstimposmean =  41;     % 41 mean of sig pos place correlations
cclustID.corrstimnegmean =  42;     % 42 mean of sig neg place correlations

cclustID.mpp =              43;     % 43 mean MPP input

cclustID.runperc =          44;     % 44 percentage of running associated cells
cclustID.runnerp =          45;     % 45 p value running
                                    % 46 p value of cell being a responder
                                    % 47 p value of cell not being a responder
cclustID.airpp =            48;     % 48 p value number of response events
                                    % 49 p value of responded stimuli
cclustID.followcell =       50;     % 50 logical if cell was in FOV of all experiments

save('cclustID','cclustID')

%% Experiment ID

if nargin >2
    cclust(1:size(caim.C,1),cclustID.expID) = expIDin;
else
    cclust(1:size(caim.C,1),cclustID.expID) = NaN;
end
%% number events and percentage of running asociated events

cclust(:,cclustID.meanf:cclustID.meanf+2) = caim.fireprop.fire;
cclust(:,cclustID.runperc) = sum(caim.S_bin(:,scn.running==1),2)./sum(caim.S_bin,2);

%% place coding strength

if isfield(scn,'cellID')  
    isplace = false(1,size(caim.C,1));
    isplace(scn.cellID(:,1)) = true;
    cclust(isplace,cclustID.plcfld) = scn.cellID(:,7);
    cclust(isplace,cclustID.plcfldp) = scn.cellID(:,3);
    cclust(scn.cellID(scn.cellID(:,2)==0,1),cclustID.plcfldp) = nan;
    cclust(~isplace,cclustID.plcfldp) = nan;
    cclust(isplace,cclustID.plcvctang) = scn.cellID(:,4);
    cclust(isplace,cclustID.plcvct) = scn.cellID(:,5);
    cclust(isplace,cclustID.plcvctp) = scn.cellID(:,6);
    cclust(~isplace,cclustID.plcvctp) = nan;
    cclust(isplace,cclustID.plclength) = scn.cellID(:,8); 
    cclust(isplace,cclustID.plcfirst) = scn.firstround;
    cclust(:,cclustID.speedcorr) = scn.speedcorr.cellID(:,1);
    cclust(:,cclustID.speedcorrp1) = scn.speedcorr.cellID(:,3);
    cclust(:,cclustID.speedcorrp2) = scn.speedcorr.cellID(:,4);
    cclust(isplace,cclustID.spatialIN) = scn.In;
    cclust(isplace,cclustID.spatialINp) = scn.cellID(:,9); 
    cclust(~isplace,cclustID.spatialINp) = nan;
else
    cclust(:,cclustID.plcfld) = NaN;
    cclust(:,cclustID.plcfldp) = NaN;
    cclust(:,cclustID.plcvctang) = NaN;
    cclust(:,cclustID.plcvct) = NaN;
    cclust(:,cclustID.plcvctp) = NaN; 
    cclust(:,cclustID.plclength) = NaN; 
    cclust(:,cclustID.plcfirst) = NaN;
    cclust(:,cclustID.speedcorr) = NaN;
    cclust(:,cclustID.speedcorrp1) = NaN;
    cclust(:,cclustID.speedcorrp2) = NaN;
    cclust(:,cclustID.spatialIN) = NaN;
    cclust(:,cclustID.spatialINp) = NaN;
end
    

%% Network event particpation

% percentage of network related/non related firing
for i = 1:size(caim.S_bin,1)
   cclust(i,cclustID.netperc) = caim.network.netprob(i)/sum(caim.S_bin(i,:));        
   cclust(i,cclustID.netnon) =  1-(caim.network.netprob(i)/sum(caim.S_bin(i,:)));
end

% normalized number of network events
cclust(:,cclustID.netprob) = caim.network.netprob./length(caim.network.netpos);

%% Correlations within networks
if isfield(caim.network,'netcorr')
coact = caim.network.coact;
cclust(:,cclustID.clust) = caim.network.cellID(:,1);
netcorr = caim.network.netcorr;
netcorr(1:size(netcorr,1)+1:end) = NaN;

for i = 1:size(coact,1)
    if size(caim.network.netraster,2)>1 && ~isempty(netcorr)
        cclust(i,cclustID.coactsum) =   nansum(coact(i,:)>1);                              % number of coactive cells
        cclust(i,cclustID.coactmean) =  nanmean(coact(i,coact(i,:)>1));                    % mean number of coactive events 
        cclust(i,cclustID.corrpossum) =  nansum(netcorr(i,:)>0);                            % number of sig pos correlated cells
        cclust(i,cclustID.corrnegsum) = nansum(netcorr(i,:)<0);                            % number of sig neg correlated cells
        cclust(i,cclustID.corrposmean) =   nanmean(netcorr(i,netcorr(i,:)>0));                % mean of sig pos correlations
        cclust(i,cclustID.corrnegmean) = nanmean(netcorr(i,netcorr(i,:)<0));                % mean of sig neg correlations        
    else
        cclust(i,cclustID.coactsum) =         NaN;
        cclust(i,cclustID.coactmean) =        NaN;
        cclust(i,cclustID.corrpossum) =       NaN;
        cclust(i,cclustID.corrposmean) =      NaN;
        cclust(i,cclustID.corrnegsum) =       NaN;
        cclust(i,cclustID.corrnegmean) =      NaN;              
    end
end
 % correlation network to place
for i = 1:size(coact,1)
    if isfield(scn,'cellID') && size(caim.network.netraster,2)>1 && ~isempty(netcorr)
        isplace(scn.cellID(scn.cellID(:,6)>.05,1)) = false;       
        cclust(i,cclustID.coactplacesum) =  sum(coact(i,isplace)>1);                % number of coactive place cells
        cclust(i,cclustID.coactplacemean) = mean(coact(i,coact(i,isplace)>1));      % mean of coactive place cell - events
        a = netcorr(i,:)>0 & isplace;
        cclust(i,cclustID.corrplacepossum) = nansum(a);              % number of sig pos correlated place cells
        cclust(i,cclustID.corrplaceposmean) = nanmean(netcorr(i,a));  % mean of sig pos place correlations
        a = netcorr(i,:)<0 & isplace;
        cclust(i,cclustID.corrplacenegsum) = nansum(a);              % number of sig neg correlated place cells
        cclust(i,cclustID.corrplacenegmean) = nanmean(netcorr(i,a));  % mean of sig neg place correlations
    else
        cclust(i,cclustID.coactplacesum) =    NaN;     
        cclust(i,cclustID.coactplacemean) =   NaN; 
        cclust(i,cclustID.corrplacepossum) =  NaN;
        cclust(i,cclustID.corrplacenegsum) =  NaN;    
        cclust(i,cclustID.corrplaceposmean) = NaN;  
        cclust(i,cclustID.corrplacenegmean) = NaN;   
    end
end
end
%% Stimulus coding strength

if isfield(scn,'airpuff')
    stim = scn.airpuff;
    isstim = false(1,size(caim.C,1));
    isstim(stim.cellID(:,1)) = true;
    % number of stimulus related events
    cclust(isstim,cclustID.nevent) = stim.cellID(:,2);
    cclust(stim.nr.cellID(:,1),cclustID.nevent) = stim.nr.cellID(:,2);
    
    % number of reponded stimuli
    cclust(isstim,cclustID.nstim) = stim.cellID(:,3);
    cclust(stim.nr.cellID(:,1),cclustID.nstim) = stim.nr.cellID(:,3);
    
    % place preference vector
    cclust(isstim,cclustID.stimang) = angle(cell2mat(stim.stimpol(:,1)));
    cclust(isstim,cclustID.stimleng) = cell2mat(stim.stimpol(:,2));
    cclust(stim.nr.cellID(:,1),cclustID.stimang) = NaN;
    cclust(stim.nr.cellID(:,1),cclustID.stimleng) = NaN;
    
    % Correlations within network events
    for i = 1:size(coact,1) 
        if size(caim.network.netraster,2)>1  && ~isempty(netcorr)
            cclust(i,cclustID.coactstimsum) = sum(coact(i,isstim)>1);                % number of coactive stim cells
            cclust(i,cclustID.coactstimmean) = mean(coact(i,coact(i,isstim)>1));      % mean of coactive place cell - events
            cclust(i,cclustID.corrstimpossum) = nansum(netcorr(i,isstim)>0);              % number of sig pos correlated place cells
            cclust(i,cclustID.corrstimnegsum) = nansum(netcorr(i,isstim)<0);              % number of sig neg correlated place cells           
            a = netcorr(i,:)>0 & isstim;
            cclust(i,cclustID.corrstimposmean) = nanmean(netcorr(i,a));  % mean of sig pos place correlations
            a = netcorr(i,:)<0 & isstim;
            cclust(i,cclustID.corrstimnegmean) = nanmean(netcorr(i,a));  % mean of sig neg place correlations
        else
            cclust(i,cclustID.coactstimsum) = NaN;
            cclust(i,cclustID.coactstimmean) = NaN;
            cclust(i,cclustID.corrstimpossum) = NaN;
            cclust(i,cclustID.corrstimnegsum) = NaN;
            cclust(i,cclustID.corrstimposmean) = NaN;
            cclust(i,cclustID.corrstimnegmean) = NaN;
        end
    end

else 
    cclust(:,cclustID.nevent) = NaN;
    cclust(:,cclustID.nstim) = NaN;
    cclust(:,cclustID.stimang) = NaN;
    cclust(:,cclustID.stimang) = NaN;
    cclust(:,cclustID.stimleng) = NaN;
end


%% Input MPP strength

if isfield(caim,'bulk')
    wind = [-.5 .5];
    Fs = 1000/scn.tsscn(end)*length(scn.tsscn);
    wind = round(wind(1)*Fs):round(wind(end)*Fs);
    for i = 1:size(caim.S,1)
        a = find(caim.S_bin(i,:));
        a = a(a>abs(wind(1)) & a < size(caim.S,2)-wind(end));
        if a
            b = zeros(length(wind),1);
            for j = 1:length(a)
                b = b+caim.bulk.trace(a(j)+wind,2);
            end
            b = b./length(a);
            cclust(i,cclustID.mpp) = max(b);
        end
    end
else
    cclust(:,cclustID.mpp) = NaN;
end

%% if there, add p values from shuffle analysis

if isfield(caim,'shuffle')
    pcell = caim.shuffle.pcell;
    cclust(:,cclustID.runnerp) = pcell(:,1);
    cclust(:,cclustID.airpp-2) = pcell(:,2); 
    cclust(:,cclustID.airpp-1) = pcell(:,3); 
    cclust(:,cclustID.airpp)   = pcell(:,4); 
    cclust(:,cclustID.airpp+1) = pcell(:,5); 
end

end