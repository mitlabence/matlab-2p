function [dist,shufdist,confint,distHist,distPerc,pInd] = PCAshuf(netraster,hasSpikesBool,score,xDim)
%%
% if sum(hasSpikesBool) <xDim; xDim = sum(hasSpikesBool);end
num = 1:xDim;

numit = 1000;
dist = cossim(netraster,num,score,hasSpikesBool);
shufdist = zeros(size(netraster,2),length(num),length(numit));
bin = .01;
x = -1:bin:1;
xx = -1-bin/2:bin:1+bin/2;

parfor i = 1:numit
    shufraster = shufshuf(netraster);
    shufdist(:,:,i) = cossim(shufraster,num,score,hasSpikesBool);
end

distcount = histcounts(dist(~isnan(dist(:))),xx,'normalization','probability');
randdist = histcounts(shufdist(~isnan(shufdist(:))),xx,'normalization','probability');

%%

a = cumsum(randdist);
b = [find(a>=.05,1) find(a>.95,1)];
confint = x([b(1) b(2)-1]);

distHist = [x;
    distcount;
    randdist];

distPerc = [sum(dist(~isnan(dist(:)))<confint(1,1))/sum(~isnan(dist(:))) sum(dist(~isnan(dist(:)))>confint(1,2))/sum(~isnan(dist(:)));
    sum(shufdist(~isnan(shufdist(:)))<confint(1,1))/sum(~isnan(shufdist(:))) sum(shufdist(~isnan(shufdist(:)))>confint(1,2))/sum(~isnan(shufdist(:)))];
%     sum(dist(~isnan(dist(:)))<confint(2,1))/length(dist(:)) sum(dist(~isnan(dist(:)))>confint(2,2))/length(dist(:));
%     sum(shufdist(:)<confint(2,1))/length(shufdist(:)) sum(shufdist(:)>confint(2,2))/length(shufdist(:))];

%% individual p values
% first frame: logical if dist is above 95% threshold of random data
% second frame: logical if the indivudual comparisson is significant
% third frame: p value for individual comparission

pInd = zeros(size(dist,1),size(dist,2),3);

A = dist;
A(A>confint(1,1)&A<confint(1,2)) = 0;
A(abs(A)>0) = 1;

B = nan(size(dist));
for i = 1:size(dist,1)
    for j = 1:size(dist,2)
        if dist(i,j)>0
            B(i,j) = sum(shufdist(i,j,:)>dist(i,j))/size(shufdist,3);
        elseif dist(i,j)<0
            B(i,j) = sum(shufdist(i,j,:)<dist(i,j))/size(shufdist,3);
        end
    end
end
C = zeros(size(B));
C(B<=.05) = 1;
pInd(:,:,1) = A;
pInd(:,:,2) = C;
pInd(:,:,3) = B;

end


function dist = cossim(netraster,num,score,hasSpikesBool)
    dist = zeros(size(netraster,2),length(num));
    for j = 1:length(num)
        k = num(j);
%         aa = zeros(size(netraster,1),1);
        aa = score(:,k);
        for i = 1:size(netraster,2)  
    %                 score = result.kern.estParams.L;
            dist(i,j) = (aa'*netraster(1:end,i))./(sqrt(sum(netraster(:,i).^2))*sqrt(sum(aa.^2)));
        end
    end   
end

function shufraster = shufshuf(netraster)
%%
    shufraster = zeros(size(netraster));
    for jj = 1:size(netraster,1)
            a = randperm(size(netraster,2));
            shufraster(jj,:) = netraster(jj,a);
    end
end
