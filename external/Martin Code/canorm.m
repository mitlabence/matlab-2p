function [S_norm,S_bin] = canorm(caim)
%%
S_norm = zeros(size(caim.S));
S_bin = zeros(size(caim.S));
for i = 1:size(caim.S,1)
%%
    b = zeros(1,length(caim.S(i,:)));
    b(caim.C(i,:)>caim.thresh(i,:)) = 1;
    b(1:end-2) = b(3:end);
    S_norm(i,:) = caim.S(i,:).*b;
    S_norm(i,:) = S_norm(i,:)/max(S_norm(i,:));
end
% caim.S_norm = caim.S_norm/max(caim.S_norm(:));
S_bin(S_norm>0) = 1 ; 

end 