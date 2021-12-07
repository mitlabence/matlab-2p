C = zeros(size(caim.C));
Cdf = zeros(size(caim.C));
Czo = zeros(size(caim.C));
for i = 1:size(caim.C,1)
    C(i,:) = caim.C(i,:);
    Cdf(i,:) = caim.C(i,:)/caim.Df(i);
%     between 0 1 
    Czo(i,:) = mat2gray(caim.C(i,:));
end


 
figure
plot(mean(Czo))


%%events quantification
Mu=mean(Czo)
sigma=std(Czo)
Threshold=Mu+5*sigma
Czobin=Czo>Threshold
plot(Czobin)


%%
figure
plot(C')
hold on
plot(mean(C),'linewidth',2)
%%
figure
plot(sum(caim.S_bin,1))
hold on
plot(sum(C),'linewidth',1)