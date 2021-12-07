function [numshuf,thresh] = ShuffleClust(caim,numit,thresh,numclustin,threshornum)
%%

network = caim.network;
netraster = network.netraster;
netnum = network.netnum;
numshuf = zeros(length(numit),1);
randr = zeros(numit,200);
% figure
% hold on

% [numclust,~] = clustsort(netraster,numclustin,thresh);
%%
if threshornum == 0
    parfor i = 1:numit
        shufraster = shufshuf(netraster);
        [numclust,~] = clustsort(shufraster,numclustin,thresh);
        numshuf(i) = numclust;
    end
    thresh = [];
else
    parfor i = 1:numit
        shufraster = shufshuf(netraster);
        randrtemp = threshfit(shufraster);
        randr(i,:) = randrtemp;
    end
    x = -1+0.005:.01:1-0.005;
    y = cumsum(sum(randr)/sum(randr(:)));
    thresh = x(find(y>0.95,1));
    numshuf = [];
end

end

function shufraster = shufshuf(netraster)
    shufraster = zeros(size(netraster));
    for j = 1:size(netraster,1)
        a = randperm(size(netraster,2));
        shufraster(j,:) = netraster(j,a);
    end
end

function randr = threshfit(netraster)
    [netcorr,~] = corr(netraster');
    a = ~isnan(diag(netcorr));
    netcorr = netcorr(a,a);
    Y = netcorr;
    Y(1:size(Y,1)+1:end) = NaN;
    randr = histcounts(Y,-1:.01:1);
end

function [newnum,clustThresh] = clustsort(netraster,numclust,thresh)
    %% cluster sorting
    [netcorr,~] = corr(netraster');    
    a = ~isnan(diag(netcorr));
    netcorr = netcorr(a,a);
    Y = netcorr;
    
%     celldist = zeros(size(netraster,1));
%     for i = 1:length(celldist)
%         for j = 1:length(celldist)
%             celldist(i,j) = (netraster(j,:)*netraster(i,:)')./(sqrt(sum(netraster(j,:)))* sqrt(sum(netraster(i,:))) );
%         end
%     end
%     netcorr = celldist;
%     Y = celldist;
    
    Z = linkage(Y,'weighted','seuclidean');
    
    % r-threshold for cluster number

    clustThresh = [];
    k = 1;
    for j = -round(numclust/2):1:round(numclust)
        C = cluster(Z,'maxclust',numclust+j);
        B = zeros(1,max(C));

        for i = 1:max(C)
            A = netcorr(C==i,C==i);
            A = tril(A,-1);
            A(A==0) = NaN;
            B(i) =  nanmean(A(:));
        end
        if nanmean(B)>thresh && k == 1
            newnum = numclust+j;
            k = 0;
        end
        clustThresh(j+round(numclust/2)+1,:) = [numclust+j,nanmean(B)];
    end
    if k == 1
        newnum = numclust+j;
    end
%     C = cluster(Z,'maxclust',numclust);
%     [~,c] = sort(C,'descend');
%     netcluster = Y(c,c);
%     imagesc(netcluster)
%     colormap(jet)

%     for i = 1:size(netraster,1)
%         for j = 1:size(netraster,1)
%             Y(i,j) = netraster(i,:)*netraster(j,:)'/(sqrt(sum(abs(netraster(i,:))))*sqrt(sum(abs(netraster(j,:)))));
%         end
%     end
    
%     plot(clustThresh(:,1),clustThresh(:,2))
end

