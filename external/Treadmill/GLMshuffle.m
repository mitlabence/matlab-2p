function [glmp,glmshuffle,compshuffle] = GLMshuffle(caim,scn,method)
if ~isfield(caim,'C') || isempty(caim.C)
    out = [];
    return
end

numit = 1000;
% method = 'pca';
mdlperf = caim.GLMout.mdlperf;
glmshuffle = zeros(10,6,2,numit);
compshuffle = zeros([size(caim.GLMout.compperf) numit]);

parfor i = 1:numit
    shufflecaim = SchuffleGLMData(caim);
    GLMout = GLMana(shufflecaim,scn,method,0);
    glmshuffle(:,:,:,i) = GLMout.mdlperf;
    compshuffle(:,:,:,:,i) = GLMout.compperf;
end
%%
glmp = sum(mdlperf>glmshuffle,4)/numit;

out.glmp = glmp;
% out.glmshuffle = glmshuffle;
% out.compshuffle = compshuffle;
out.glmshuffleMean  = mean(glmshuffle,4);
out.glmshuffleStd   = std(glmshuffle,0,4);
out.compshuffleMean = mean(compshuffle,5);
out.compshuffleStd  = std(compshuffle,0,5);

% out.glmp = caim.GLMshuffle.glmp;
% out.glmshuffleMean = mean(caim.GLMshuffle.glmshuffle,4);
% out.glmshuffleStd = std(caim.GLMshuffle.glmshuffle,0,4);  
% out.compshuffleMean = mean(caim.GLMshuffle.compshuffle,5);
% out.compshuffleStd = std(caim.GLMshuffle.compshuffle,0,5); 

end

function DataOut = SchuffleGLMData(DataIn)

caim = DataIn;
%%
C = caim.C;
% shuffle everything
%     for j = 1:size(S_bin,1)
%         a = randperm(size(S_bin,2));
%         caim.S_bin(j,:) = S_bin(j,a);
%         caim.C(j,:) = caim.C(j,a);
%     end

%     shift traces randomly with respect to each other
for j = 1:size(C,1)
    a = randperm(size(C,2),1);
    C(j,:) = C(j,[a:end 1:a-1]);
end
%%
caim.C = C;
DataOut = caim;

end