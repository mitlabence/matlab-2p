function caim = nikonCutToLabView(belt_struct,caim)
%NIKONCUTTOLABVIEW If Nikon stopped recording after LabView, this
%function cuts the Nikon recording to match LabView.
%   Copied from Martin's readcaim.m function.
% Input:
%   belt_struct: belt data in a struct (e.g. correctBeltInplace.m)
%   caim: Nikon opened as caim struct (CNMF)
if ~isfield(belt_struct, 'tsscn')
   error("Error in matchNikonToLabView: belt_struct parameter has no property 'tsscn'. Possible cause: incorrect type of belt_struct. Use correctBeltInplace()!");
end
tsscn = belt_struct.tsscn;

if isfield(caim,'Y') && length(tsscn)<size(caim.Y,2)
    disp('Shortening following to match scanner time frame: Y, C, S, f, thresh, S_norm, S_bin');
    disp(['reason: tsscn field (' length(tsscn) ') is shorter than Y field (' size(caim.Y,2) ').']);
    disp('Possibly dangerous: calling loadArray of CNMF object to change Y (data) parameter. It should, however, update Yr and the new dimensions.');
    caim.loadArray(caim.Y(:,1:length(tsscn)));
    %caim.Y = caim.Y(:,1:length(tsscn)); %change back to this in case
    %loadArray does not work. But this is a very bad way to change data!
    %Need to change d1, d2, Yr...
    caim.C = caim.C(:,1:length(tsscn));
    caim.S = caim.S(:,1:length(tsscn));
    caim.f = caim.f(:,1:length(tsscn));
    caim.thresh = caim.thresh(:,1:length(tsscn));
    caim.S_norm = caim.S_norm(:,1:length(tsscn));
    caim.S_bin = caim.S_bin(:,1:length(tsscn));
else
    disp('CMNF object data was not shortened.');
end
end

