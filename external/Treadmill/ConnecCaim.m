function [belt,caim] = ConnecCaim(pathname,files,samedate)

if exist('CAIM','var');clear CAIM;end
if exist('BELT','var');clear BELT;end
    
for i = 1:length(samedate)  
    filename = files(samedate(i)).name;
    filename = filename(1:end-6);
    [belt,caim] = readcaim(pathname,filename); 
    if isfield(caim,'C');[caim.S_bin,caim.S_norm] = normCa(caim);end
    CAIM(i) = caim;
    BELT(i) = belt;  
end

if length(samedate) == 1
    return
else
    disp([num2str(length(samedate)) ' data sets are merged to one experiment'])
end

%% Concatenate belt Data

belt.tsscn = BELT(1).tsscn;
belt.round = BELT(1).round;
belt.speed = BELT(1).speed;
belt.distance = BELT(1).distance;
belt.distancePR = BELT(1).distancePR;
belt.reflect = BELT(1).reflect;
belt.licking = BELT(1).licking;
belt.stripes = BELT(1).stripes;
belt.stripesPR = BELT(1).stripesPR;
belt.time = BELT(1).time;
belt.timePR = BELT(1).timePR;
belt.reward = BELT(1).reward;
belt.airpuff = BELT(1).airpuff;
belt.soundl = BELT(1).soundl;
belt.soundr = BELT(1).soundr;

if isfield(BELT,'odor1')
    belt.odor1 = BELT(1).odor1;
    belt.odor2 = BELT(1).odor2;
    belt.odor3 = BELT(1).odor3;
    belt.odor4 = BELT(1).odor4;
    belt.odor5 = BELT(1).odor5;
end

if isfield(BELT,'pupil')
    belt.pupil = BELT(1).pupil;
end

tcut = zeros(4,2);
for i = 2:length(BELT)

    % find matching position on the belt between two experiments
    if BELT(i).distancePR(1) < belt.distancePR(end)
        lastrnd = belt.distancePR(belt.round == belt.round(end));
%           figure;plot(lastrnd);hold on;plot(BELT(i).distancePR(BELT(i).round==0));hold off;
    else
        lastrnd = belt.distancePR(belt.round == belt.round(end) | belt.round == belt.round(end)-1);
%             figure;plot(lastrnd);hold on;plot(BELT(i).distancePR(BELT(i).round==0));hold off;
    end
    tcut(i,1) = find(abs(lastrnd - BELT(i).distancePR(1)) == min(abs(lastrnd - BELT(i).distancePR(1))),1);
    tcut(i,1) = length(lastrnd)-tcut(i,1);
    tcut(i,2) = find(abs(belt.tsscn-belt.time(end-tcut(i,1))) == min(abs(belt.tsscn-belt.time(end-tcut(i,1)))));
    tcut(i,2) = length(belt.tsscn)-tcut(i,2);

    belt.tsscn = belt.tsscn(1:end-tcut(i,2)); 
    belt.tsscn = [belt.tsscn; BELT(i).tsscn + max(belt.tsscn)+mean(diff(belt.tsscn))];
    belt.round = belt.round(1:end-tcut(i,1)); 
    belt.round = [belt.round; BELT(i).round+max(belt.round)];

    belt.speed = belt.speed(1:end-tcut(i,1)); 
    belt.speed = [belt.speed; BELT(i).speed];
    belt.distance = belt.distance(1:end-tcut(i,1)); 
    belt.distance = [belt.distance; BELT(i).distance+belt.distance(end)];        
    belt.distancePR = belt.distancePR(1:end-tcut(i,1)); 
    belt.distancePR = [belt.distancePR; BELT(i).distancePR];

    belt.reflect = belt.reflect(1:end-tcut(i,1)); 
    belt.reflect = [belt.reflect; BELT(i).reflect];
    belt.licking = belt.licking(1:end-tcut(i,1)); 
    belt.licking = [belt.licking; BELT(i).licking];
    belt.stripes = belt.stripes(1:end-tcut(i,1)); 
    belt.stripes = [belt.stripes; BELT(i).stripes+belt.stripes(end)];
    belt.stripesPR = belt.stripesPR(1:end-tcut(i,1)); 
    belt.stripesPR = [belt.stripesPR; BELT(i).stripesPR];
    belt.time = belt.time(1:end-tcut(i,1)); 
    belt.time = [belt.time; BELT(i).time+belt.time(end)+mean(diff(belt.time))];        
    belt.timePR = belt.timePR(1:end-tcut(i,1)); 
    belt.timePR = [belt.timePR; BELT(i).timePR];
    belt.reward = belt.reward(1:end-tcut(i,1)); 
    belt.reward = [belt.reward; BELT(i).reward];
    belt.airpuff = belt.airpuff(1:end-tcut(i,1)); 
    belt.airpuff = [belt.airpuff; BELT(i).airpuff];
    belt.soundl = belt.soundl(1:end-tcut(i,1)); 
    belt.soundl = [belt.soundl; BELT(i).soundl];
    belt.soundr = belt.soundr(1:end-tcut(i,1)); 
    belt.soundr = [belt.soundr; BELT(i).soundr];

    if isfield(BELT,'odor1')
        belt.odor1 = belt.odor1(1:end-tcut(i,1));
        belt.odor1 = [belt.odor1; BELT(i).odor1];
        belt.odor2 = belt.odor2(1:end-tcut(i,1));             
        belt.odor2 = [belt.odor2; BELT(i).odor2];
        belt.odor3 = belt.odor3(1:end-tcut(i,1)); 
        belt.odor3 = [belt.odor3; BELT(i).odor3]; 
        belt.odor4 = belt.odor4(1:end-tcut(i,1)); 
        belt.odor4 = [belt.odor4; BELT(i).odor4];
        belt.odor5 = belt.odor5(1:end-tcut(i,1)); 
        belt.odor5 = [belt.odor5; BELT(i).odor5];
    end

    if isfield(BELT,'pupil')
        belt.pupil = belt.pupil(1:end-tcut(i,1)); 
        belt.pupil = [belt.pupil; BELT(i).pupil];
    end
end

%%
if isfield(caim,'Y')
    Cm = cell(length(CAIM),1);
    numcell = zeros(1,length(samedate));
    % Correlate alle reference picture in case of drift during the experiment
    [optimizer, metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 0.001;
    optimizer.Epsilon = 1.5e-5;
    optimizer.GrowthFactor = 1.01;
    optimizer.MaximumIterations = 300;
    ref = CAIM(1).Cn;
    refpic = 1;

    for j = 1:length(CAIM)

        CAIM(j).tform = imregtform(CAIM(j).Cn,ref, 'rigid', optimizer, metric);
        TF(j).tform = CAIM(j).tform;
        ref = CAIM(j).Cn;  

        i = j-1;
        while i>refpic      
            CAIM(j).tform.T = CAIM(j).tform.T*TF(i).tform.T;
            i = i-1;
        end


        CAIM(j).Cncor = imwarp(CAIM(j).Cn,CAIM(j).tform,'OutputView',imref2d(size(ref)));
        tform = imregtform(CAIM(j).Cncor,CAIM(refpic).Cn, 'affine', optimizer, metric);

        CAIM(j).tform.T = CAIM(j).tform.T*tform.T;
        CAIM(j).Cncor = imwarp(CAIM(j).Cn,CAIM(j).tform,'OutputView',imref2d(size(ref)));                


        B = CAIM(j).A;
        d1 = size(CAIM(j).Cn,1);
        d2 = size(CAIM(j).Cn,2);
        BB = zeros(size(B));
        c = zeros(1,size(B,2));
        for i = 1:size(B,2)
            b = full(reshape(B(:,i),[d1 d2]));
            b = imwarp(b,CAIM(j).tform,'OutputView',imref2d(size(ref)));
            if find(b,1);c(i) = 1;  end
            BB(:,i) = reshape(b,[d1*d2 1]);
        end
        BB = sparse(BB(:,c==1));
        CAIM(j).A = BB;
        CAIM(j).cm = com(BB,d1,d2);  
    end


      
    for i = 1:length(CAIM)
        A = CAIM(i).A;
        d1 = CAIM(i).options.d1;
        d2 = CAIM(i).options.d2;
        cm = com(A,d1,d2);
        Cm{i} = cm; 
        numcell(i) = size(cm,1);         
    end
   
    indcell = 1:length(Cm);
end
%% Finding Cells with close distance between recordings and high x-correllation 

if isfield(caim,'Y')
    [samecell,corrcell] = findsamecell(CAIM,Cm,indcell);

    for i = 2:length(indcell)-1
        indcell1 = indcell(i:end);
        indcell1(indcell1>indcell(i-1)) = indcell1(indcell1>indcell(1))-(i-1);
        Cm1 = Cm(indcell(i:end));

        [samecell1,corrcell1] = findsamecell(CAIM(indcell(i:end)),Cm1,1:size(Cm1,1));

        unused = ones(1,size(samecell1,2));
        unused(samecell(indcell(i),samecell(indcell(i),:)~=0)) = 0;
        samecell2 = zeros(size(samecell,1),size(samecell1,2));
        corrcell2 = samecell2;
        samecell2(indcell(i:end),:) = samecell1;
        corrcell2(indcell(i:end),:) = corrcell1;
        samecell2 = samecell2(:,unused==1);
        corrcell2 = corrcell2(:,unused==1);
        samecell = [samecell samecell2];
        corrcell = [corrcell corrcell2];
    end


    Cm1 = Cm(indcell(end));
    unused = ones(1,size(Cm1{1},1));
    unused(samecell(indcell(end),samecell(indcell(end),:)~=0)) = 0;
    samecell2 = zeros(size(samecell,1),length(find(unused)));
    corrcell2 = samecell2;
    samecell2(indcell(end),:) = find(unused);
    samecell = [samecell samecell2];
    corrcell = [corrcell corrcell2];
end
%% Concatenate traces of same cells or fill passive time series
t = zeros(1,length(CAIM)+1);

for i = 1:length(CAIM)
    t(i) = t(i)-tcut(i,2);%evtl +1...
    t(1) = 0;
    if isfield(caim,'Y')
        t(i+1) = t(i)+size(CAIM(i).Y,2); 
    elseif isfield(caim,'bulk')
         t(i+1) = t(i)+size(CAIM(i).bulk.traceMEC,2); 
    end
end

if isfield(caim,'Y')
    caim.Y = zeros(size(samecell,2),t(end));
    caim.A = zeros(size(CAIM(1).A,1),size(samecell,2));
    caim.C = caim.Y;
    caim.S = caim.Y;
    caim.b = CAIM(1).b;
    caim.f = zeros(1,t(end));   %!!!
    caim.Cn = CAIM(1).Cn;
    caim.Df = zeros(size(samecell,2),1);
    caim.options = CAIM(1).options;
    caim.cID = ones(size(samecell,2),1);
    caim.thresh = caim.Y;   
    caim.S_norm = caim.Y;
    caim.S_bin = caim.Y;

    for i = 1:size(samecell,2)
        iii = 0;    
        for ii = 1:size(samecell,1)
            if samecell(ii,i) > 0
                caim.Y(i,t(ii)+1:t(ii+1)) = CAIM(ii).Y(samecell(ii,i),1:t(ii+1)-t(ii));
                caim.C(i,t(ii)+1:t(ii+1)) = CAIM(ii).C(samecell(ii,i),1:t(ii+1)-t(ii));%/CAIM(ii).Df(samecell(ii,i));
                caim.S(i,t(ii)+1:t(ii+1)) = CAIM(ii).S(samecell(ii,i),1:t(ii+1)-t(ii));
                caim.thresh(i,t(ii)+1:t(ii+1)) = CAIM(ii).thresh(samecell(ii,i),1:t(ii+1)-t(ii));
                caim.S_norm(i,t(ii)+1:t(ii+1)) = CAIM(ii).S_norm(samecell(ii,i),1:t(ii+1)-t(ii));
                caim.S_bin(i,t(ii)+1:t(ii+1)) = CAIM(ii).S_bin(samecell(ii,i),1:t(ii+1)-t(ii));
                if iii == 0
                    caim.A(:,i) = CAIM(ii).A(:,samecell(ii,i));
                    caim.Df(i) = CAIM(ii).Df(samecell(ii,i));
                    iii = 1;
                end
            end
        end       
    end
end

if isfield(CAIM,'bulk')
   if isfield(CAIM(1).bulk,'cropCA1')
       caim.bulk.cropCA1 = CAIM(1).bulk.cropCA1;
       caim.bulk.threshCA1 = CAIM(1).bulk.threshCA1;
       caim.bulk.traceCA1 = zeros(2,t(end)); 
       caim.bulk.templateCA1 = CAIM(1).bulk.templateCA1;
   end
   if isfield(CAIM(1).bulk,'cropMEC')
       caim.bulk.cropMEC = CAIM(1).bulk.cropMEC;
       caim.bulk.threshMEC = CAIM(1).bulk.threshMEC;
       caim.bulk.traceMEC = zeros(2,t(end)); 
       caim.bulk.templateMEC = CAIM(1).bulk.templateMEC; 
   end
   if isfield(CAIM(1).bulk,'cropGC')
       caim.bulk.cropGC = CAIM(1).bulk.cropGC;
       caim.bulk.threshGC = CAIM(1).bulk.threshGC;
       caim.bulk.traceGC = zeros(2,t(end)); 
       caim.bulk.templateGC = CAIM(1).bulk.templateGC; 
   end
end

for i = 1:size(CAIM,2)
    if isfield(CAIM,'f')
        caim.f(t(i)+1:t(i+1)) = CAIM(i).f(1:t(i+1)-t(i));
    end
    if isfield(CAIM,'bulk')
        if isfield(CAIM(1).bulk,'cropCA1')
        caim.bulk.traceCA1(:,t(i)+1:t(i+1)) = CAIM(i).bulk.traceCA1(:,1:t(i+1)-t(i));
        end
        if isfield(CAIM(1).bulk,'cropMEC')
        caim.bulk.traceMEC(:,t(i)+1:t(i+1)) = CAIM(i).bulk.traceMEC(:,1:t(i+1)-t(i));  
        end
        if isfield(CAIM(1).bulk,'cropGC')
        caim.bulk.traceGC(:,t(i)+1:t(i+1)) = CAIM(i).bulk.traceGC(:,1:t(i+1)-t(i));
        end
    end
end


end

function [samecell,corrcell] = findsamecell(CAIM,Cm,indcell)

cm = Cm{indcell(1)};
d1 = CAIM(1).options.d1;
d2 = CAIM(1).options.d2;
maxdist = CAIM(1).options.gSig*2;
samecell = zeros(length(CAIM),size(cm,1));
samecell(indcell(1),:) = 1:length(cm);
corrcell = zeros(length(CAIM),size(cm,1));


thresh = .8;
for kk = 2:length(CAIM)
    k = indcell(kk);
    cm2 = Cm{k};   
    for i = 1:length(cm)
        distcell = zeros(1,length(cm2));            
        for jj = 1:length(cm2)
            distcell(jj) = norm(cm2(jj,:)-cm(i,:));
        end
        samec = find(distcell < maxdist);
        if isempty(samec)
            samecell(k,i) = 0;
        else
            corrc = zeros(1,length(samec));
            for ii = 1:length(samec)
                a = full(reshape(CAIM(indcell(1)).A(:,i),d1,d2));
                b = full(reshape(CAIM(k).A(:,samec(ii)),d1,d2));
                corrc(ii) = max(max(xcorr2(a,b)));
            end
            if max(corrc(:))>=thresh
                [corrcell(k,i),a] = max(corrc(:));
                b = find(samecell(k,:)==samec(a));
                if isempty(b)
                    samecell(k,i) =samec(a);
                elseif corrcell(k,b) < corrcell(k,i)                        
                    samecell(k,i) =samec(a);
                    samecell(k,b) = 0;
                end
            else
                samecell(k,i) =0;
            end                
        end        
    end    
end

%%
nr = size(samecell,2);
sx = min([CAIM(1).options.sx,floor(d1/2),floor(d2/2)]);
int_x = zeros(nr,2*sx);
int_y = zeros(nr,2*sx);
start = 11;
stop = start+10;

h = figure('position',[680 49 1016 948]);
for i = start:stop     
    for j = 1:length(CAIM)
        if samecell(j,i)>0
        subplot(stop-start+1,length(CAIM),j+(i-start)*length(CAIM))
        a = reshape(CAIM(j).A(:,samecell(j,i)),d1,d2);
        int_x(i,:) = round(Cm{j}(samecell(j,i),1)) + (-(sx-1):sx);
        if int_x(i,1)<1
            int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
        end
        if int_x(i,end)>d1
            int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
        end
        int_y(i,:) = round(Cm{j}(samecell(j,i),2)) + (-(sx-1):sx);
        if int_y(i,1)<1
            int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
        end
        if int_y(i,end)>d2
            int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
        end
        a = a(int_x(i,:),int_y(i,:));
        fa = imagesc(int_x(i,:),int_y(i,:),a); axis square;   
        end
    end
end

end
    
function [S_bin,S_norm] = normCa(caim)
%%
thresh = zeros(size(caim.C));
sigmamult = 6;
bins = 50;
win = 6000;        % Decrease window size for descrete adaptation of threshold (smooth threshold scroll down)
numwin = ceil(size(caim.C,2)/win);
counts = zeros(numwin,bins);
centers = zeros(numwin,bins);
threshtemp = zeros(numwin,1);
gaussEqn = 'a*exp(-((x-b)/(sqrt(2)*c))^2)';

for i = 1:size(caim.S,1)
    
    %%
    for j = 1:numwin
        %%       
        if j < numwin
            tempwin = 1+(j-1)*win:j*win;
        else
            tempwin = 1+(j-1)*win:size(caim.S,2);
        end
        a = caim.S(i,tempwin);
        [counts(j,:),hh] = histcounts(a(a~=0),bins);
        centers(j,:) = hh(1:end-1);
    end
    
    %%
    for j = 1:numwin
        x = centers(j,:);
        y = counts(j,:);
        countmax = find(y==max(y),1);
        
        f1 = fit(x',y',gaussEqn,...
            'Start', [y(countmax) x(countmax) 10],...
            'upper',[1000 0 1000],...
            'lower',[0 0 0]);
        threshtemp(j) = f1.b+sigmamult*abs(f1.c);
        
%         bar(x,y)
%         hold on
%         plot(x,f1(x))
%         hold off
    end

    %%
    for j = 1:numwin
        if j < numwin
            tempwin = 1+(j-1)*win:j*win;
        else
            tempwin = 1+(j-1)*win:size(caim.S,2);
        end
        thresh(i,tempwin) = threshtemp(j);
    end
    
   
end

%% Smooth threshold adaptation following photobleaching (exponatial fit)

% a = sum(caim.S,1)';
% a = a-min(a);
% a = a/max(a);
% d = belt.tsscn/1000;
% fs = length(d)/d(end);  
% fc = [1]; % cutoff frequency
% Wn= fc/(fs/2);
% n = 2;
% [b,bb] = butter(n,Wn,'low');
% aa = filtfilt(b,bb,a);
% d = d(1:length(aa));
% 
% ExpEqn = 'a*exp(-((x-b)/c))+d';
% startPoints = [.2 0 1000 min(aa)];
% upperBounds = [1 1 10000 .4];
% lowerBounds =[.05 -1 500 0];
% f1 = fit(d,aa,ExpEqn,'Start', startPoints,'upper',upperBounds,'lower',lowerBounds,'Exclude',d<100 & d>1000);
% aa = f1(d);
% 
% thresh=thresh.*(1+20*aa');


%%
S_norm = zeros(size(caim.S));
S_bin = zeros(size(caim.S));
for i = 1:size(caim.S,1)

    S_bin(i,caim.S(i,:)>thresh(i,:)) = 1;
%     plot(b);hold on
%     b(1:end-2) = b(3:end);
%     plot(b);hold off
    
    S_norm(i,:) = caim.S(i,:).*S_bin(i,:);
    S_norm(i,:) = S_norm(i,:)/max(S_norm(i,:));
end
% 
% S_bin(S_norm>0) = 1; 

%% sparsen inter spike intervall

S_sparse = S_bin;
for j = 1:size(S_bin,1)   
    a = find(S_bin(j,:));
    S_sparse(j,a(diff(a)==1)+1) = 0;
end

S_bin = S_sparse;
S_norm = S_bin.*S_norm;

%% output
caim.thresh = thresh;
caim.S_bin = S_bin;
caim.S_norm = S_norm;
end
