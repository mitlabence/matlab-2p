%  Load pooled Data
load('/media/2Photon/Nicola/BigFatSummary.mat')


%% Correlate Experiments
% all sessions
sessions = 1:8;
refpic = [1 1 1 1 1 1 1 1 1 1];
ID = { '177' '184' '235' '255' '339zone1' '339zone2' '342' '349new'  };

cclust = [];
samecell = [];
plcfield = [];
coact = [];
if exist('CAIMcorr','var'); clear CAIMcorr;end
%%
for i = [ 3]%1:size(CAIM,2)
    %%
    useexp = false(1,size(CAIM,1));
    for j = sessions%1:length(useexp)
        if ~isempty(CAIM(j,i).behave) %ischar(ID{i})
            useexp(j) = 1;
        end
    end
    %%
    [A,cclusttemp,samecelltemp,plcfieldtemp,coacttemp,CAIMcorr(useexp==1,i)] = FollowCells(CAIM(useexp==1,i),ID{i},refpic(i));
    cclust(end+1:end+size(cclusttemp,1),:,useexp==1) = cclusttemp;
    samecell(useexp==1,end+1:end+size(samecelltemp,2)) = samecelltemp;
    plcfield(end+1:end+size(cclusttemp,1),:,useexp==1)  = plcfieldtemp;
%     coact(end+1:end+size(cclusttemp,1),:,useexp==1) = coacttemp;
end

%%
save(['/media/2Photon/Nicola/BigFatCorrData.mat'],'CAIMcorr','ID','cclust','samecell','plcfield')






  






