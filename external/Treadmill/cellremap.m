function remapOUT = cellremap(CAIMcorr,mouse,expID)

experiment = {'Base1','Base2','Base3','Base4','Base5','Cues1','Cues2','Cues3','Air1','Air2','Air3','AirFix','Retr'};

load('cclustID.mat','cclustID');

%%
emptyCAIM = false(1,length(expID));
for i = 1:length(expID)
    if isempty(CAIMcorr(expID(i),mouse).A)
        emptyCAIM(i) = true;
    end
end
expID(emptyCAIM) = [];
cclust = cell(length(expID),1);
netraster = cell(length(expID),1);
netnum = zeros(length(expID)+1,1);
isplace = cell(length(expID),1);
isstim = cell(length(expID),1);
isboth = cell(length(expID),1);
isthere = cell(length(expID),1);

for i = 1:length(expID)
    trial = expID(i);
    netraster{i} = CAIMcorr(trial,mouse).network.netraster;
    netnum(i+1) = size(netraster{i},2)+netnum(i);
    cclusttemp = CAIMcorr(trial,mouse).cclust;
    
%     active field sorting
    netraster{i}(CAIMcorr(trial,mouse).inField==0,:) = 0;
    cclusttemp(CAIMcorr(trial,mouse).inField==0,:) = NaN;

    isplace{i} = (cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & ~(cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
    isstim{i}  = ~(cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & (cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
    isboth{i}  = (cclusttemp(:,cclustID.plcvct)>0 & cclusttemp(:,cclustID.plcvctp)<=.05) & (cclusttemp(:,cclustID.nstim)>2 & cclusttemp(:,cclustID.airpp)<=.05);
    isthere{i} = ~isnan(cclusttemp(:,1));
    cclust{i} = cclusttemp;
end

%% fill up individual indicators so that they all have same length

for i = 1:length(expID-1)
    netraster{i}(size(netraster{end},1),:) = 0;
    isplace{i}(size(netraster{end},1)) = 0;
    isstim{i}(size(netraster{end},1)) = 0;
    isboth{i}(size(netraster{end},1)) = 0;
    isthere{i}(size(netraster{end},1)) = 0;
    cclust{i}(size(netraster{end},1),1) = NaN;
end

%%
cclusttemp = cell2mat(cclust(1));
for i = 2:length(cclust)
    cclusttemp = cat(3,cclusttemp,cell2mat(cclust(i)));
end
cclust = cclusttemp;

% transfer the individual cells to one matrix with aditional indicator in
% first coloumn
isplace = cell2mat(isplace');
isplace(:,2:end+1) = isplace; 
isplace(:,1) = sum(isplace(:,2:end),2);
isplace(isplace>1) = 1;

isstim = cell2mat(isstim');
isstim(:,2:end+1) = isstim; 
isstim(:,1) = sum(isstim(:,2:end),2);
isstim(isstim>1) = 1;

isboth = cell2mat(isboth');
isboth(:,2:end+1) = isboth; 
isboth(:,1) = sum(isboth(:,2:end),2);
isboth(isboth>1) = 1;

isstim(isboth(:,1)==1,1)=0;
% isstim(isplace(:,1)==1,1)=0;
isplace(isboth(:,1)==1,1)=0;
% isplace(isstim(:,1)==1,1)=0;

isthere = cell2mat(isthere');
isthere(:,2:end+1) = isthere;

%% properties of remapping cells

cellID = isplace(:,1)>0 | isboth(:,1)>0;
cclusttemp = cclust(cellID,:,:);
cclusttemp = permute(cclusttemp,[1 3 2]);
angles = cclusttemp(:,:,cclustID.plcvctang);
angles(cclusttemp(:,:,cclustID.plcvct)<.1) = NaN;
shiftangles = angles;
shiftangles(angles<0) = shiftangles(angles<0)+2*pi;
remap = angles;
numfields = zeros(1,size(remap,1));

for i = 1:size(angles,1)
    a = find(~isnan(angles(i,:)));
    numfields(i) = 1;
    remap(i,a(1)) = numfields(i);
    for j = 1:length(a)-1
        k = j; % Check all fields until first one is found
        while k >0
            if angles(i,a(j+1))>=angles(i,a(k))-pi/4 && angles(i,a(j+1))<=angles(i,a(k))+pi/4
               remap(i,a(j+1)) =remap(i,a(k));
               k = 0;
            elseif shiftangles(i,a(j+1))>=shiftangles(i,a(k))-pi/4 && shiftangles(i,a(j+1))<=shiftangles(i,a(k))+pi/4
               remap(i,a(j+1)) =remap(i,a(k));
               k = 0;
            elseif k == 1
               numfields(i) = numfields(i) +1;
               remap(i,a(j+1)) = numfields(i);
               k = 0;
            else 
               k = k-1;
            end
        end
    end
end

%%


field = false(size(remap));
fieldold =  false(size(remap));
fieldnew = false(size(remap));

for i = 1:size(remap,1)
    %%
    if ~isnan(remap(i,1))
        field(i,1) = true;
    end
    for j = 2:size(remap,2)
        if ~isnan(remap(i,j)) && sum(~isnan(remap(i,1:j-1)))>0
            if sum(remap(i,j) == remap(i,1:j-1))>0
                fieldold(i,j) = true;
            else
                fieldnew(i,j) = true;
            end
        elseif ~isnan(remap(i,j))
            field(i,j) = true;
        end
    end
end

remapOUT = NaN(3,length(emptyCAIM));
remapOUT(:,~emptyCAIM) = [sum(field,1); sum(fieldold,1); sum(fieldnew,1)];

end