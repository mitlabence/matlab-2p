function [A,cclust,samecell,plcfield,coact,CAIM] = FollowCells(CAIM,mouseID,refpic)
%%

[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.001;
optimizer.Epsilon = 1.5e-5;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;

%% check if somatic data is present
if refpic == 0 || (isempty(CAIM(1).Cn) && isempty(CAIM(1).A))
    A = [];
    cclust = [];
    samecell = [];
    plcfield = [];
    coact = [];
    for j = 1:length(CAIM)
        a = zeros(100);
        a(100) = 1;
        CAIM(j).tform  =  [];%imregtform(a,a, 'rigid', optimizer, metric);
        CAIM(j).Cncor = [];
        CAIM(j).B = [];
        CAIM(j).cm = [];
        CAIM(j).ActField = [];
        CAIM(j).inField = [];
    end
    return
end


%%

ref = CAIM(1).Cn;
Cm = cell(length(CAIM),1);
numcell = zeros(1,length(CAIM));

for j = 1:length(CAIM)
    if ~isempty(CAIM(j).Cn)
        
        CAIM(j).tform = imregtform(CAIM(j).Cn,ref, 'rigid', optimizer, metric);
        TF(j).tform = CAIM(j).tform;
        ref = CAIM(j).Cn;  

        i = j-1;
        while i>refpic      
            CAIM(j).tform.T = CAIM(j).tform.T*TF(i).tform.T;
            i = i-1;
        end


        CAIM(j).Cncor = imwarp(CAIM(j).Cn,CAIM(j).tform,'OutputView',imref2d(size(ref)));
        tform = imregtform(CAIM(j).Cncor,CAIM(refpic).Cn, 'affine', optimizer, metric);

        CAIM(j).tform.T = CAIM(j).tform.T*tform.T;
        CAIM(j).Cncor = imwarp(CAIM(j).Cn,CAIM(j).tform,'OutputView',imref2d(size(ref)));                


        B = CAIM(j).A;
        d1 = size(CAIM(j).Cn,1);
        d2 = size(CAIM(j).Cn,2);
        BB = zeros(size(B));
        c = zeros(1,size(B,2));
        for i = 1:size(B,2)
            b = full(reshape(B(:,i),[d1 d2]));
            b = imwarp(b,CAIM(j).tform,'OutputView',imref2d(size(ref)));
            if find(b,1);c(i) = 1;  end
            BB(:,i) = reshape(b,[d1*d2 1]);
        end
        BB = sparse(BB(:,c==1));
        CAIM(j).B = BB;
        CAIM(j).cclust = CAIM(j).cclust(c==1,:);
        CAIM(j).plcfield = CAIM(j).plcfield(c==1,:);
        CAIM(j).cm = com(BB,d1,d2);
        Cm{j} = CAIM(j).cm;
        numcell(j) = size(Cm{j},1);   
    else
        CAIM(j).B = [];
        CAIM(j).cclust = [];
        CAIM(j).plcfield = [];
        Cm{j} = [];
        numcell(j) = 0; 
    end
end

indcell = 1:length(Cm);


%%


cm = Cm{indcell(1)};
A = CAIM(indcell(1)).B;
d1 = size(CAIM(1).Cn,1);
d2 = size(CAIM(1).Cn,2);
maxdist = 5;%CAIM(indcell(1)).options.gSig*2;
samecell = zeros(length(CAIM),size(cm,1));
samecell(indcell(1),:) = 1:length(cm);
corrcell = zeros(length(CAIM),size(cm,1));

cclust = CAIM(1).cclust;
if size(CAIM(1).plcfield,2)<100
    CAIM(1).plcfield = NaN(size(CAIM(1).A,2),100);
end
plcfield = CAIM(1).plcfield(:,1:100);
coact = CAIM(1).network.coact;

for kk = 2:length(CAIM)
    k = indcell(kk);
    cm2 = Cm{k};   
    for i = 1:length(cm)
        distcell = zeros(1,length(cm2));            
        for jj = 1:length(cm2)
            distcell(jj) = norm(cm2(jj,:)-cm(i,:));
        end
        samec = find(distcell < maxdist);
        if isempty(samec)
            samecell(k,i) = 0;
        else
            corrc = zeros(1,length(samec));
            for ii = 1:length(samec)
                a = full(reshape(A(:,i),d1,d2));
                b = full(reshape(CAIM(k).B(:,samec(ii)),size(CAIM(k).Cn)));
                corrc(ii) = max(max(xcorr2(a,b)));
            end
            thresh =.5;
            if max(corrc(:))>=thresh
                [corrcell(k,i),a] = max(corrc(:));
                b = find(samecell(k,:)==samec(a));
                if isempty(b)
                    samecell(k,i) =samec(a);
                elseif corrcell(k,b) < corrcell(k,i)                        
                    samecell(k,i) =samec(a);
                    samecell(k,b) = 0;
                end
            else
                samecell(k,i) =0;
            end                
        end      
    end
    
    %% look for new cells in 2nd data set
    
    newcell = 1:size(cm2,1);
    newcell(samecell(k,(samecell(k,:)~=0))) = [];  
    A = [A CAIM(k).B(:,newcell)];
    cm = [cm; cm2(newcell,:)];      
    samecell(k,end+1:end+length(newcell)) = newcell;
    
    % spatial components
    CAIM(k).B(:,samecell(k,:)~=0) = CAIM(k).B(:,samecell(k,samecell(k,:)~=0));
    CAIM(k).B(:,samecell(k,:)==0) = 0;
    
    % center of mass
    Cm{k}(samecell(k,:)~=0,:) = Cm{k}(samecell(k,samecell(k,:)~=0),:); 
%     CAIM(k).cm(samecell(k,:)~=0,:) =Cm{k}(samecell(k,samecell(k,:)~=0),:);
%     CAIM(k).cm(samecell(k,:)==0,:) = NaN;
    CAIM(k).cm = cm;

    % ccluster
    cclust(samecell(k,:)~=0,:,k) = CAIM(k).cclust(samecell(k,samecell(k,:)~=0),:);
    CAIM(k).cclust(samecell(k,:)~=0,:) = CAIM(k).cclust(samecell(k,samecell(k,:)~=0),:);
    CAIM(k).cclust(samecell(k,:)==0,:) = NaN;
    
    % placefields
    if size(CAIM(k).plcfield,2)<100
        CAIM(k).plcfield(samecell(k,:)~=0,:) = NaN(size(CAIM(k).A,2),100);
    end
    plcfield(samecell(k,:)~=0,:,k) = CAIM(k).plcfield(samecell(k,samecell(k,:)~=0),1:100);
    CAIM(k).plcfield(samecell(k,:)~=0,1:100) = CAIM(k).plcfield(samecell(k,samecell(k,:)~=0),1:100);
    
    % networks
    CAIM(k).network.netraster(samecell(k,:)~=0,:) = CAIM(k).network.netraster(samecell(k,samecell(k,:)~=0),:);
    CAIM(k).network.netraster(samecell(k,:)==0,:) = NaN;
    CAIM(k).network.netprob(samecell(k,:)~=0) = CAIM(k).network.netprob(samecell(k,samecell(k,:)~=0));
    for i = 1:length(CAIM(k).network.netID)        
        CAIM(k).network.netID(i) = {find(CAIM(k).network.netraster(:,i))};
    end
    CAIM(k).network.coact(samecell(k,:)~=0,samecell(k,:)~=0) = CAIM(k).network.coact(samecell(k,samecell(k,:)~=0),samecell(k,samecell(k,:)~=0));
    CAIM(k).network.coact(samecell(k,:)==0,samecell(k,:)==0) = NaN; 
    
    if ~isempty(CAIM(k).network.netID) && size(CAIM(k).network.netraster,2)>2 && isfield(CAIM(k).network,'netcorr')
        CAIM(k).network.netcorr(samecell(k,:)~=0,samecell(k,:)~=0) = CAIM(k).network.netcorr(samecell(k,samecell(k,:)~=0),samecell(k,samecell(k,:)~=0));
        CAIM(k).network.netcorr(samecell(k,:)==0,samecell(k,:)==0) = NaN; 
    else
        CAIM(k).network.netcorr = NaN(size(samecell));
    end

end

%% Active window
    
c = ones(size(CAIM(1).Cn));
for k = 1:size(CAIM)
    a = ones(size(CAIM(k).Cn));
    b = imwarp(a,CAIM(k).tform,'OutputView',imref2d(size(a)));
    c = c.*b;
end

% Define for individual experiment if cells are in the FOV
for k = 1:size(CAIM)
    inField = zeros(1,size(Cm{k},1));
    CAIM(k).ActField = c;
    for j = 1: size(Cm{k},1)
        if round(Cm{k}(j,1)) > 0 && round(Cm{k}(j,2)) > 0 && round(Cm{k}(j,1)) < size(c,1) && round(Cm{k}(j,2)) < size(c,2)
            inField(j) = c(round(Cm{k}(j,1)),round(Cm{k}(j,2)));
        end
    end
    CAIM(k).inField = inField;
    CAIM(k).cclust(:,end+1) = inField;
end

% Define for the comlete set if cells are in the FOV
inField = zeros(size(cm,1),1,size(CAIM,1));
for j = 1: size(cm,1)
    if round(cm(j,1)) > 0 && round(cm(j,2)) > 0 && round(cm(j,1)) < size(c,1) && round(cm(j,2)) < size(c,2)
        inField(j,1,:) = c(round(cm(j,1)),round(cm(j,2)));
    end
end
cclust(:,end+1,:) = inField;
%% Plot results

tempdir  = '/media/2Photon/Nicola/temp/';
figure('color',[0 0 0],...
    'position',[500 50 1.5*[420 594]],...
    'renderer','painters',...
    'visible','off')

j = 1;
k = 1;
for i = 1:size(CAIM,1)
    if j > 4*3
%         exportfig([tempdir mouseID '-' num2str(k)])
        printpdf([tempdir mouseID '-' num2str(k)])
        k = k+1;
        j = 1;
        figure('color',[0 0 0],...
            'position',[500 50 1.5*[420 594]],...
            'visible','off')
    end
    subplot(4,3,j)
    imshowpair(CAIM(i).Cn,CAIM(i).Cn,'Scaling','joint')
    
    subplot(4,3,j+1)
    imshowpair(CAIM(refpic).Cn,CAIM(i).Cncor,'Scaling','joint')
        
    subplot(4,3,j+2)
        
%     imshowpair(reshape(sum(CAIM(i).B,2),size(CAIM(i).Cn)),reshape(sum(CAIM(1).B,2),size(CAIM(1).Cn)),'Scaling','joint')
    imshowpair(full(reshape(sum(CAIM(i).B,2),size(CAIM(i).Cn))),full(reshape(sum(CAIM(refpic).B,2),size(CAIM(1).Cn))),'Scaling','joint')
    j = j+3; 
end
% exportfig([tempdir mouseID '-' num2str(k)])
printpdf([tempdir mouseID '-' num2str(k)])
%%
close all
if isfile(['/media/2Photon/Nicola/CellTracing/' mouseID  '-CellCorr.pdf'])
    delete(['/media/2Photon/Nicola/CellTracing/' mouseID  '-CellCorr.pdf'])
end
files = dir([tempdir '*.pdf']);
olddir = cd;
cd(tempdir)
append_pdfs(['/media/2Photon/Nicola/CellTracing/' mouseID  '-CellCorr.pdf'],files.name)
delete(files.name)
cd(olddir)


end

