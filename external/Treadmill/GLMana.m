function GLMout = GLMana(caim,scn,method,runshuffle)
if ~isfield(caim,'C')
    GLMout = [];
    return
end
addpath(genpath('/media/2Photon/Matlab code/Martin Source Code/gpfa_v0203'))


C = caim.C; C = C./caim.Df;
% C = caim.S_bin;C(caim.S_bin==0)= 0;
C(:,scn.running==0) = 0;
% kernSD = 5;
% C = smoother(C, kernSD,1);
% Include only place cells
% C = C(scn.cellID(scn.cellID(:,6)<=.05,1),:);

meanC       = mean(C,2);
C_centered  = C-repmat(meanC,1,size(C,2));
stdC        = std(C,[],2);
C_normed    = C_centered./repmat(stdC+.001,1,size(C,2));
C = C_normed;

    
if nargin>2 && strcmp(method,'ica')
    %%
    if nargin<4 || runshuffle ~= 0;disp('GLM using ICA components');end
    xDim = 30;
    yOut = smoother(C, kernSD,1);
    b = rica(yOut',xDim,'IterationLimit',1000);
    % b = rica(a,q,'IterationLimit',1000);
    a = yOut'*b.TransformWeights;
    score = b.TransformWeights;
    hasSpikesBool = true(size(C,1),1);
else
    if nargin<4 || runshuffle ~= 0;disp('GLM using PCA components');end
    %%  
    % Extract neural trajectories
%     [U,S,V]=svd(C_normed);
%     a = C_normed'*U;
    a = pca(C);
end

%% linear model fit

speed = scn.speed*100;
speed(speed<0) = 0;
speed = smooth(speed,30);

running = scn.running;
distance = scn.distance;

alpha = linspace(-pi,pi,150);
distance(distance>1500) = 1500;
theta = zeros(length(distance),1);
sinth = zeros(length(distance),1);
costh = zeros(length(distance),1);
for i = 1:length(scn.distance)
    if scn.distance(i) >0
        theta(i) = alpha(ceil(distance(i)/10));
        sinth(i) = sin(theta(i));
        costh(i) = cos(theta(i));
    else
        theta(i) = 0;
        sinth(i) = 0;
        costh(i) = 1;
    end
end

movement = [sinth costh speed running theta];
brdln = round(4/5*size(a,1));
ypred = zeros(size(movement));
test = zeros(size(movement));
for i = 1:4
    if i<4
        % gaussian process regression        
%         mdl1(i).mdl = fitrgp(a(1:brdln,:),movement(1:brdln,i));       
%         mdl1(i+size(movement,2)).mdl = fitrgp(C(:,1:brdln)',movement(1:brdln,i));

        % binary regression tree
%         mdl1(i).mdl = fitrtree(a(1:brdln,:),movement(1:brdln,i));       
%         mdl1(i+size(movement,2)).mdl = fitrtree(C(:,1:brdln)',movement(1:brdln,i));

        % generalized linear model
        mdl1(i).mdl = fitglm(a(1:brdln,:),movement(1:brdln,i));       
        mdl1(i+size(movement,2)).mdl = fitglm(C(:,1:brdln)',movement(1:brdln,i));
        
        ypred(:,i) = predict(mdl1(i).mdl,a);
        test(:,i) = ypred(:,i)-movement(:,i);   
        ypred(:,i+size(movement,2)) = predict(mdl1(i+size(movement,2)).mdl,C');
        test(:,i+size(movement,2)) = ypred(:,i+size(movement,2))-movement(:,i); 
    else        
        mdl1(i).mdl = fitclinear(a(1:brdln,:),movement(1:brdln,i));
        mdl1(i+size(movement,2)).mdl = fitclinear(C(:,1:brdln)',movement(1:brdln,i));
        ypred(:,i) = predict(mdl1(i).mdl,a);
        test(:,i) = abs(ypred(:,i)-movement(:,i));   
        ypred(:,i+size(movement,2)) = predict(mdl1(i+size(movement,2)).mdl,C');
        test(:,i+size(movement,2)) = abs(ypred(:,i+size(movement,2))-movement(:,i)); 
    end
      
end

% reconstruct position from cos & sin

ypred(:,5) = angle(1j*(ypred(:,1))+ypred(:,2));
ypred(:,10) = angle(1j*(ypred(:,6))+ypred(:,7));
test(:,5) = angle(pi+ypred(:,2)+1j*(ypred(:,1))-(costh+1j*sinth));
test(:,10) = angle(pi+ypred(:,7)+1j*(ypred(:,6))-(costh+1j*sinth));

% intervalls for training and testing: complete, running and resting 
int = true(1,size(a,1)); int(1,1:brdln) = 0;
int(2,:) = scn.running == 1; int(2,1:brdln) = 0;
int(3,:) = scn.running == 0; int(3,1:brdln) = 0;
int(4,:) = true(1,size(a,1)); int(4,brdln+1:size(a,1)) = 0;
int(5,:) = scn.running == 1; int(5,brdln+1:size(a,1)) = 0;
int(6,:) = scn.running == 0; int(6,brdln+1:size(a,1)) = 0;

% read out the mean and std of deviation from actual movement for two types of input, 
% for all three parameters 
% in all six time intervalls

mdlperf = zeros(size(ypred,2),size(int,1),2);
for i = 1:size(ypred,2)
    for j = 1:size(int,1)
        mdlperf(i,j,1:2) = [mean(test(int(j,:),i)) std(test(int(j,:),i))];
    end
end

GLMout.movement = movement;
GLMout.pred = ypred;
GLMout.test = test;
GLMout.mdlperf = mdlperf;
%% Succesive increasing number of components
% runIdx = 2;
% method = 'pca';
% xDim = size(C,1)-1;
% binWidth = 1; 
% kernSD = 5;
% 
% % Extract neural trajectories
% result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
%                     'kernSDList', kernSD,'binWidth',binWidth);
% if isempty(result);PCAout = [];return;end
% [estParams, seqTrain] = postprocess(result, 'kernSD', kernSD);%%
% score = result.kern.estParams.L;
% hasSpikesBool = result.hasSpikesBool;
% if sum(hasSpikesBool)<xDim
%     xDim = sum(hasSpikesBool);
% end
% %    score = estParams.Corth;
% a = [];
% for i = 1:size(seqTrain,2)
%     a = [a seqTrain(i).xorth];
% end
% clear b
% 
% for i = 1:size(a,1)
%     b(i,:) = interp1(1:size(a,2),a(i,:),1:size(a,2)/size(caim.C,2):size(a,2),'spline');
% end
% 
% b(:,end+1:size(caim.C,2)) = 0;
% a = b';

%%
maxnum = 100;
num = linspace(log(3),log(maxnum),15);
num = [round(exp(num(1:end)))];
if size(a,2)<maxnum;num = num(num<size(a,2));end

mdlperf = zeros(size(ypred,2),size(int,1),2,length(num));
ypred = zeros(size(movement));
test = zeros(size(movement));

clear mdl
for jj = 1:length(num)
    numint = 1:num(jj);
    for i = 1:4
        if i<4
            mdl1(i).mdl = fitlm(a(1:brdln,numint),movement(1:brdln,i));
            mdl1(i+size(movement,2)).mdl = fitglm(C(numint,1:brdln)',movement(1:brdln,i));   
            ypred(:,i) = predict(mdl1(i).mdl,a(:,numint));
            test(:,i) = ypred(:,i)-movement(:,i);
            ypred(:,i+size(movement,2)) = predict(mdl1(i+size(movement,2)).mdl,C(numint,:)');
            test(:,i+size(movement,2)) = ypred(:,i+size(movement,2))-movement(:,i); 
        else        
            mdl1(i).mdl = fitclinear(a(1:brdln,numint),movement(1:brdln,i));
            mdl1(i+size(movement,2)).mdl = fitclinear(C(numint,1:brdln)',movement(1:brdln,i));
            ypred(:,i) = predict(mdl1(i).mdl,a(:,numint));
            test(:,i) = abs(ypred(:,i)-movement(:,i));
            ypred(:,i+size(movement,2)) = predict(mdl1(i+size(movement,2)).mdl,C(numint,:)');
            test(:,i+size(movement,2)) = abs(ypred(:,i+size(movement,2))-movement(:,i));       
        end       
    end

    % reconstruct position from cos & sin

    ypred(:,5) = angle(1j*(ypred(:,1))+ypred(:,2));
    test(:,5) = angle(pi+ypred(:,2)+1j*(ypred(:,1))-(costh+1j*sinth));
    ypred(:,10) = angle(1j*(ypred(:,6))+ypred(:,7));
    test(:,10) = angle(pi+ypred(:,7)+1j*(ypred(:,6))-(costh+1j*sinth));
    for i = 1:size(ypred,2)
        for j = 1:size(int,1)
            mdlperf(i,j,1:2,jj) = [mean(test(int(j,:),i)) std(test(int(j,:),i))];
        end
    end
end

GLMout.numcomp = num;
GLMout.compperf = mdlperf;

%% shuffle analysis
if isfield(caim,'GLMout') && isfield(caim.GLMout,'glmp')
    GLMout.glmp = caim.GLMout.glmp;
    GLMout.glmshuffleMean  = caim.GLMout.glmshuffleMean;
    GLMout.glmshuffleStd   = caim.GLMout.glmshuffleStd;
    GLMout.compshuffleMean = caim.GLMout.compshuffleMean;
    GLMout.compshuffleStd  = caim.GLMout.compshuffleStd;    
end

if nargin<4
    runshuffle = false;
end

if runshuffle == 1
    disp('run GLM shuffle')
    caim.GLMout = GLMout;
    [glmp,glmshuffle,compshuffle] = GLMshuffle(caim,scn,method);
    
    GLMout.glmp = glmp;
    GLMout.glmshuffleMean  = mean(glmshuffle,4);
    GLMout.glmshuffleStd   = std(glmshuffle,0,4);
    GLMout.compshuffleMean = mean(compshuffle,5);
    GLMout.compshuffleStd  = std(compshuffle,0,5);    
else


end