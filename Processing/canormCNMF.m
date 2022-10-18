function [S_norm,S_bin] = canormCNMF(caim, thresh)
%CANORMCNMF copy of canorm from Martin's code, suitable for CNMF objects.
%It does not change the CNMF object (caim). Where to get thresh from?
% example: [cID, thresh] = sort_components(CNM); (CNM is a CNMF object.)
% sort_components is Martin's code. cID is used in
% CNM = divcellsCNMF(CNM, cID, thresh);
% but this renders CNM an invalid CNMF object!

disp('Function canormCNMF has been called. This function is an undocumented copy of canorm. It should be refactored!');
S_norm = zeros(size(caim.S));
S_bin = zeros(size(caim.S));
%%
for i = 1:size(caim.S,1) %go through each neuron. S size is #neurons x #frames
    b = zeros(1,length(caim.S(i,:))); %vector of size 1 x #frames
    b(caim.C(i,:)>thresh(i,:)) = 1; %binary vector, temporal component of neuron is greater than provided threshold
    b(1:end-2) = b(3:end); %shift 2 to the left? why?
    S_norm(i,:) = caim.S(i,:).*b;
    S_norm(i,:) = S_norm(i,:)/max(S_norm(i,:));
end
% caim.S_norm = caim.S_norm/max(caim.S_norm(:));
S_bin(S_norm>0) = 1 ; %could just return b instead, too...

end 
