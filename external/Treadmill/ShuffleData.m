

numit = 1000;
S_bin = caim.S_bin;
C = caim.C;
Df = zeros(size(C));
Df(:,1) = caim.Df;
netevents(caim,scn);
networkshuffle = zeros(numit,3);
fireshuffle = zeros(numit,3);
%% inter spike intervall

S_sparse = S_bin;
for j = 1:size(S_bin,1)
    
    a = find(S_bin(j,:));
    S_sparse(j,a(diff(a)==1)+1) = 0;
%     plot(S_bin(j,:));hold on;plot(S_sparse(j,:)+1);hold off
end

network = netevents(cat(3,S_sparse,S_sparse),scn);
disp(size(network.netpos,2))
S_bin = S_sparse;
%%

parfor i = 1:numit
    
    S_temp = zeros(size(S_bin));
    C_temp = zeros(size(S_bin));
    % shuffle everything
    for j = 1:size(S_bin,1)
        a = randperm(size(S_bin,2));
        S_temp(j,:) = S_bin(j,a);
        C_temp(j,:) = C(j,a);
    end
    
    % shift traces randomly with respect to each other
%     for j = 1:size(S_bin,1)
%         a = randperm(size(S_bin,2),1);
%         S_temp(j,:) = S_bin(j,[a:end 1:a-1]);
%         C_temp(j,:) = C(j,[a:end 1:a-1]);
%     end
    
%%

    network = netevents(cat(3,S_temp,C_temp),scn);
    if ~isempty(network.netfreq)
        networkshuffle(i,:) = network.netfreq(1,:);
    end
    
    fireprop = firefreq(cat(3,S_temp,C_temp,Df),scn);
    fireshuffle(i,:) = fireprop.meanfire(2,:);
    
%     [belt,scn] = stimcor(belt,caimshuffle,scn);
%     scn.cclust = cellcluster(caimshuffle,scn);
    

end

%%
figure
subplot(3,3,1)
histogram(networkshuffle(:,1))

subplot(3,3,2)
histogram(networkshuffle(:,2))

subplot(3,3,3)
histogram(networkshuffle(:,3))

% subplot(3,3,4)
% histogram(fireshuffle(:,1))

subplot(3,3,5)
histogram(fireshuffle(:,2))

subplot(3,3,6)
histogram(fireshuffle(:,3))
% imagesc(caimshuffle.network.netraster);
% plot(sum(caim.S_bin))