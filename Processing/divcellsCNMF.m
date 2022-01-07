function [caim, cID, thresh] = divcellsCNMF(caim,cID, thresh)  
%DIVCELLSCNMF this function performs the same functionality as divcells.m
%(written by Martin(?)), with one difference:
%   1. as the input object caim is a CNMF object, it cannot have the 'raw'
%   field. In divcells, this 'raw' field is a copy of the caim object. This
%   is better achieved by saving the result of this function in a new
%   object. (The input object caim is not modified, as far as I know, a
%   local copy is made within the function. I.e.
%   caim2 = divcellsCNMF(caim, cID); leaves caim unchanged, and creates the
%   caim2 object that is the result of this function.
%
% Input:
%   caim: CNMF object
%   cID: cell ID object (see sort_components.m)
%   thresh: thresh object (see sort_components.m)
% Output:
%   caim: CNMF object
%   cID: cID cut
%   thresh: thresh cut

%TODO: decode sort_components.m, finish documentation on cID here. What
%does this function do?

%% Check if input is a CNMF object
if isfield(caim,'raw')
    error('divcellsCNMF: caim object has field "raw", i.e. is not a CMNF object!');
end

%% Perform the function
caim.A = caim.A(:,cID==1);
caim.C  = caim.C(cID==1,:);
caim.Df = caim.Df(cID==1);
caim.S  = caim.S(cID==1,:);
caim.Y  = caim.Y(cID==1,:);
cID = cID(cID==1,:);
thresh = thresh(cID==1,:);
end

