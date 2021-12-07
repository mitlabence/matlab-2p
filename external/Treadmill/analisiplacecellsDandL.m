% load ID file where the property IDs are found
 load('cclustID.mat')
% load the cell cluster variable from DAta set (Or from pooling CAIM)
cclust = scn.cclust;

% Do some logic to identify losoncy cells

 isplacelos = cclust(:,cclustID.plcvctp)<.05 & cclust(:,cclustID.plcvct)>0;
% isplacelos = cclust(:,cclustID.plcvct)>0;
% Read out a property of cells e.g. vectorlength
vectorlength = cclust(isplacelos,cclustID.plcvct);




% Do the same for dombeck: 1. is there an entry? Is it significant?
% isplacedom = cclust(:,cclustID.plcfldp)<.05 & cclust(:,cclustID.plcfld)>0;

% Read out a property of cells e.g. Gaussian lenght
% gausslength = cclust(isplacedom,cclustID.plclength);

 m=mean(vectorlength)

%% all gaussian
load('cclustID.mat')
cclust = scn.cclust;


isplacedom = cclust(:,cclustID.plcfld)>0;
gausslength = cclust(isplacedom,cclustID.plclength);

%%
num = [1:8];
MeanLos = nan(length(num),size(CAIM,2));
MeanDom = nan(length(num),size(CAIM,2));
for i = 1:length(num)
    
    
    for j = 1:size(CAIM,2)
        if ~isempty(CAIM(num(i),j).cclust)
            cclust = CAIM(num(i),j).cclust;
            % Do some logic to identify losoncy cells
            isplacelos = cclust(:,cclustID.plcvctp)<.05 & cclust(:,cclustID.plcvct)>0;

            % Read out a property of cells e.g. vectorlength
            vectorlength = cclust(isplacelos,cclustID.plcvct);
            MeanLos(i,j) = mean(vectorlength);


            % Do the same for dombeck: 1. is there an entry? Is it significant?
            isplacedom = cclust(:,cclustID.plcfldp)<.05 & cclust(:,cclustID.plcfld)>0;

            % Read out a property of cells e.g. Gaussian lenght
            gausslength = cclust(isplacedom,cclustID.plclength);
            MeanDom(i,j) = mean(gausslength);
        end
    end
end
