function scn = placecor(caim,scn,numit)
%%
% Input:
% caim          Calcium imaging DAta from NNMF
% scn           Behaviour data converted to scanner time frame
% numit         Number of iterations for shuffle analysis
%
% Output:
%
% scn.spaccode  number of runnin related events per cell in bin belt space
% scn.spcsort   same as spacecode with only active cells
% scn.plcfield  smoothed place activity (first row), with dombeck criteria (second row)
% scn.spcpol    losonscy place preference: 1. row: angle of pol-angle, 2. row length of angle, 3. alpha rho  bin
% scn.shufpol   
% scn.shufield  
% scn.cellID    1. CellID, 
%               2. center of place field, 
%               3. pValue placefield, 
%               4. angle of place vector,
%               5. length of place vector, 
%               6. pValue place vector
%               7. Place field position with gauss fit
%               8. Place field width with gauss fit (1 sigma)
% scn.spcon     1. t of all running related events
%               2. t of all events of cells with significant place vector
%               3. t of all events of cells with significant place field                
%               4. t of all passes of significant place fields
% scn.speedcorr 1. CellID
%               2. bins for speedbinning (min:max, 20 bins)
%               3. Traces of Amp vs bin
%               4. Onsets of significant events at certain speed
%               5. Std of traces


disp(['The mouse ran ' num2str(max(scn.rounds)) ' rounds'])

%% mean speed values in space bins

spcbin = 10; % size of bin in mm
numbin = ceil(max(scn.distance)/spcbin);
strtrnd = 1;%double(scn.distance(1)>scn.distance(end));
speedbin = zeros(max(scn.rounds),numbin);
waitbin = zeros(max(scn.rounds),numbin);
for i = strtrnd:max(scn.rounds)
%     win =  scn.rounds==i;
    runningtemp =  scn.running(scn.rounds==i)==1;
    disttemp = scn.distance(scn.rounds==i);
    speedtemp = scn.speed(scn.rounds==i);  
    for j = 1:numbin
        int = find(disttemp>(j-1)*spcbin+1 & disttemp<j*spcbin);  
        if ~isempty(int)
            temp = speedtemp(int);
            temp = temp(temp>0);
            speedbin(i,j) = mean(temp);
            waitbin(i,j) = sum(runningtemp(int)==0);
        end
    end
end


scn.speedbin = speedbin;
scn.waitbin = waitbin;

%% Correlate activity ca response in bin on belt
minround = 5;
minev = 3;

if ~isfield(caim,'C') || max(scn.rounds)<minround
    scn.spcsort = [];
    scn.spcpol = [];
    scn.plcfield = [];
    scn.speedcorr = [];
%             scn.centerfield = [];
    disp('No sufficient number of rounds during session')
    return
end

%% speed coding

[trace,traceErr,onsets,c,p,bin] = SpeedCorr(caim,scn);

if nargin<3
    numit = 1000;
end
cshuffle = zeros(length(c),numit);
pshuffle = zeros(length(p),numit);
if ispc
    for i = 1:numit
        caim2 = ShuffleData(caim);
        [~,~,~,cshuffle(:,i),pshuffle(:,i),~] = SpeedCorr(caim2,scn);
    end
else
    parfor i = 1:numit
        caim2 = ShuffleData(caim);
        [~,~,~,cshuffle(:,i),pshuffle(:,i),~] = SpeedCorr(caim2,scn);
    end
end

for i = 1:length(c)
    p(i,2) = sum(pshuffle(i,:)<p(i,1))/numit;
    p(i,3) = sum(abs(cshuffle(i,:))>abs(c(i)))/numit;
end

speedcorr.cellID = [c p];
speedcorr.bin    = bin;
speedcorr.trace  = trace;
speedcorr.onsets = onsets;
speedcorr.error  = traceErr;
scn.speedcorr    = speedcorr;

%%
[spccode,spcsort,numbin,cellID] = CaOnPlace(caim,scn,spcbin,minev);

scn.spaccode = spccode;                     
scn.spcsort = spcsort;


if isempty(scn.spcsort)
    scn.spcpol = [];
    scn.plcfield = [];
%             scn.centerfield = [];
    disp('No place related firing')
    return
end

%% Dombeck 2010/2015 placefield

plfdlg = 150; % placefieldlength in mm to look for responses
plfdlg = ceil(plfdlg/2/spcbin); %in bins and half

[plcfield,centerfield,~] = dombeck(spcsort,plfdlg);
scn.plcfield = plcfield;           
% scn.centerfield = centerfield;
% Center of placefield
cellID(:,2) = round(centerfield*spcbin);

%% Danielson/Losonczy 2016 place preference vector

% Normalizing vector with the numbers of running related recording frames per bin
o = zeros(numbin,1); % Normalizing vector with the numbers of running related recording frames per bin
disttemp = scn.distance(scn.running==1);
for j = 1:numbin
    o(j) = length(find(disttemp>(j-1)*spcbin+1 & disttemp<j*spcbin));           
end
rnd = find(diff(scn.rounds));

[spcpol,theta,rho,firstround] = losonczy(spcsort,numbin,rnd,o);  

% size of place field (gauss fit)
fldlength = fieldlength(scn,spcbin);
spcpol(:,4) = fldlength(:,1);

% angle of place preference
cellID(:,4) = angle(cell2mat(spcpol(:,1)));
% length of place coding vector
cellID(:,5) = cell2mat(spcpol(:,2));
% position of strongest place field
cellID(:,7:8) = cell2mat(fldlength(:,2));

%% Skaggs 1993 spatial information

disttemp = scn.distance(scn.running==1);
% number of bins
bins = [2 4 5 10 20 25 100];
spcbins = floor(max(scn.distance)./bins);
spc = cell(1,length(bins));
% Variable for spatial information
In = zeros(size(cellID,1),size(bins,2));
for i = 1:length(bins)
    % firing rate per bin
    [~,spctemp,~,~] = CaOnPlace(caim,scn,spcbins(i),minev,cellID(:,1));
    % compute skaggs information
    In(:,i) = skaggs(spctemp,bins(i),spcbins(i),disttemp);
    % save spctemp for later use in shuffle analysis
    spc{i} = spctemp;
end
%%

scn.spcpol = spcpol;
scn.theta = theta;
scn.rho= rho;
scn.firstround = firstround;

%% Shuffle analysis

shufpol = zeros(size(spcsort,1),2,numit);
shufield = zeros(size(spcsort,1),2,numit);

if ispc
    for i = 1:numit      
        shufspc = zeros(size(spcsort));
% %          shuffle acitivity within each round 
%         for j = 1:size(spcsort,1) %loop over cells
%             for k = 1:size(spcsort,3) %loop over rounds
%                 shufspc(j,:,k) = spcsort(j,randperm(size(spcsort,2)),k);
%             end
%         end
% 
% %         shuffle activity over all rounds
        for j = 1:size(spcsort,1) %loop over cells
            a = reshape(spcsort(j,:,:),[1, size(spcsort,2)*size(spcsort,3)]);
            a = a(randperm(size(a,2)));
            shufspc(j,:,:) = reshape(a,[1,size(spcsort,2),size(spcsort,3)]);                        
        end

        % shuffle activity over complete measurement
%         caim2 = ShuffleData(caim);
%         [~,shufspc,~,~] = CaOnPlace(caim2,scn,spcbin,minev,cellID(:,1));

        % read out vector length for shuffled data
        [spcpol,~,~] = losonczy(shufspc,numbin,rnd,o);    
        temp = cell2mat(spcpol(:,1)); 
        temp(:,2) = cell2mat(spcpol(:,2));
        shufpol(:,:,i) = temp;
        
        % read out if PC according to arguments was found with shuffled data
        [~,centerfield] = dombeck(shufspc,plfdlg);                        
        centerfield(:,2) = centerfield;
        centerfield(centerfield(:,1)>0,2) = 1;    % logical if PF is present    
        shufield(:,:,i) = centerfield;
        
    end
else
    parfor i = 1:numit
        shufspc = zeros(size(spcsort));
%         for j = 1:size(spcsort,1)
%             for k = 1:size(spcsort,3)
%                 shufspc(j,:,k) = spcsort(j,randperm(size(spcsort,2)),k);
%             end
%         end
% 
        for j = 1:size(spcsort,1) %loop over cells
            a = reshape(spcsort(j,:,:),[1, size(spcsort,2)*size(spcsort,3)]);
            a = a(randperm(size(a,2)));
            shufspc(j,:,:) = reshape(a,[1,size(spcsort,2),size(spcsort,3)]);                        
        end

%         caim2 = ShuffleData(caim);
%         [~,shufspc,~,~] = CaOnPlace(caim2,scn,spcbin,minev,cellID(:,1));
        
        % read out cevotr length for shuffled data
        [spcpol,~,~] = losonczy(shufspc,numbin,rnd,o);    
        temp = cell2mat(spcpol(:,1)); 
        temp(:,2) = cell2mat(spcpol(:,2));
        shufpol(:,:,i) = temp;
        
        % read out if PC according to arguments was found with shuffled data
        [~,centerfield] = dombeck(shufspc,plfdlg);                        
        centerfield(:,2) = centerfield;
        centerfield(centerfield(:,1)>0,2) = 1;    % logical if PF is present    
        shufield(:,:,i) = centerfield;
        
   end
end
 
% p-value place field
cellID(:,3) = sum(shufield(:,2,:),3)./numit;
% p-value place preference
cellID(:,6) = 0;
for i = 1:size(cellID,1)
    cellID(i,6) = sum(shufpol(i,2,:)<inf & shufpol(i,2,:)>scn.spcpol{i,2})./numit;
end
%% shuffle and p-value spatial information
InShuf = zeros(size(In,1),size(In,2),numit);
if ispc 
    for j = 1:numit
        Intemp = zeros(size(In));
        for i = 1:length(bins)
            spctemp = spc{i};
            for k = 1:size(spcsort,1) %loop over cells
                a = reshape(spctemp(k,:,:),[1, size(spctemp,2)*size(spctemp,3)]);
                a = a(randperm(size(a,2)));
                spctemp(k,:,:) = reshape(a,[1,size(spctemp,2),size(spctemp,3)]); 
            end
            Intemp(:,i) = skaggs(spctemp,bins(i),spcbins(i),disttemp);
        end
        InShuf(:,:,j) = Intemp;
    end
else
    parfor j = 1:numit
        Intemp = zeros(size(In));
        for i = 1:length(bins)
            spctemp = spc{i};
            for k = 1:size(spcsort,1) %loop over cells
                a = reshape(spctemp(k,:,:),[1, size(spctemp,2)*size(spctemp,3)]);
                a = a(randperm(size(a,2)));
                spctemp(k,:,:) = reshape(a,[1,size(spctemp,2),size(spctemp,3)]); 
            end
            Intemp(:,i) = skaggs(spctemp,bins(i),spcbins(i),disttemp);
        end
        InShuf(:,:,j) = Intemp;
    end
end

% Go thru every cell individually and calculate max info and p-value
Inout = zeros(size(In,1),1);
for i = 1:size(cellID,1)
    % Skaggs information baseline corrected with the mean of the shuffled analysis
    In(i,:)  = In(i,:) - mean(InShuf(i,:,:),3);
    % Shuffled distribution is baseline correctd, too
    InShuf(i,:,:) = InShuf(i,:,:)-mean(InShuf(i,:,:),3);
    % Check the binning that gives maximum information 
    [~,a] = max(In(i,:));
    % p-Value is calculated as the fraction of info being smaller then
    % shuffled info
    cellID(i,9) = sum(In(i,a)<InShuf(i,a,:))/numit;
    Inout(i) = In(i,a);
end
%%
scn.shufpol = shufpol;
scn.shufield = shufield;
scn.In = Inout;
scn.InShuf = InShuf;
scn.cellID = cellID;

%% Timepoints of place related events and placefields

% all place related events
spcon = sum(caim.S_bin(cellID(:,1),:),1);
% all significant place vecotr cell events
spcon(2,:) = sum(caim.S_bin(cellID(cellID(:,6)<=0.05,1),:),1);
%  all place field cell events
spcon(3,:) = sum(caim.S_bin(cellID(cellID(:,2)>0 & cellID(:,3)<=0.05),:),1);           
spcon(spcon>1) = 1;
spcon(:,scn.running==0) = 0;


%% timepoints of (significant) placefield passages
plcfield = cellID(cellID(:,2)>0 & cellID(:,3)<=0.05,2);
spcon(4,:) = 0;
for i = 1:length(rnd) % loop over all rounds
    for j = 1:length(plcfield) % loop over all place fields                   
        k = find(abs(scn.distance(scn.rounds==(i-1))-plcfield(j))==min(abs(scn.distance(scn.rounds==(i-1))-plcfield(j))));                    
        if length(k) == 1
            kk = find(scn.rounds==i-1,1);
            spcon(4,k+kk) = 1;
        end
    end
end

scn.spcon = spcon;
         
% plotting
% if isfield(caim,'Df') && ~isempty(scn.spcsort) 
%     polarplots(scn)
% end



end

function [spccode,spcsort,numbin,cellID] = CaOnPlace(caim,scn,spcbin,minev,cellID)
    numbin = ceil(max(scn.distance)/spcbin);
    spccode = zeros(length(caim.Df),numbin,max(scn.rounds));
    strtrnd = double(scn.distance(1)>scn.distance(end));


    %% 
    for j = 1 : size(caim.S_bin,1)   
        for i = strtrnd:max(scn.rounds)
            disttemp = scn.distance(scn.running==1 & scn.rounds==i);
        %                catemp = caim.S_act(j,scn.running==1 & scn.rounds==i)/caim.Df_act(j);
            catemp = caim.S_bin(j,scn.running==1 & scn.rounds==i);
            if i<max(scn.rounds)
                for k = 1:numbin
                    int = find(disttemp>(k-1)*spcbin+1 & disttemp<k*spcbin);  
                    if ~isempty(int)
                        spccode(j,k,i+1) = sum(catemp(int));
                    end
                end
            elseif strtrnd == 0
                for k = 1:floor(scn.distance(1)/spcbin)
                    int = find(disttemp>(k-1)*spcbin+1 & disttemp<k*spcbin);  
                    if ~isempty(int)
                        spccode(j,k,i+1) = mean(catemp(int));
                    end
                end
            end
        end
    end
    
    if nargin<5
        spccode(spccode>1) = 1;
        a = spccode;
        a = sum(a,2);
        a(a>1) = 1;
        a = sum(a,3);
        cellID = find(a>minev);
    end
    
     spcsort = spccode(cellID,:,:); 
end


function [spcpol,theta,rho,firstround] = losonczy(spcsort,numbin,rnd,o)
%%

spcpol = cell(size(spcsort,1),3);
Alpha = linspace(0,2*pi,numbin*1);
theta = linspace(0,2*(length(rnd)+1)*pi,numbin*(1+length(rnd)));
rho = (1:length(theta))*(length(rnd)+1)/length(theta);   


for i = 1:size(spcsort,1)   % loop over cells                     
    c = [];
    d = [];
    dd = [];
    for j = 1:size(spcsort,3)  % loop over rounds
          b = find(spcsort(i,:,j));
          for jj = 1:length(b)
              c = [c; Alpha(b(jj)),rho(numbin*(j-1)+b(jj)),j, spcsort(i,b(jj),j)];                        
%                           d = [d;o(b(jj))*exp(Alpha(b(jj))*1j)];
              d = [d;exp(Alpha(b(jj))*1j)/o(b(jj))];
%                           d = [d;1/o(b(jj))];
              dd = [dd; o(b(jj))];
          end  
    end
% dd/mean(dd)
% d = sum(d)*sum(dd);

    spcpol{i,1} = sum(dd)*sum(d)/length(d)^2; % place coding vector
    spcpol{i,2} = abs(sum(dd)*sum(d)/length(d)^2); % place coding vector length
    spcpol{i,3} = c; % alpha % rho % round % count in bin
%                 [length(d) max(dd) mean(dd) sum(dd)/length(d)] 
    
end
plfdlg = round(150/1500*size(spcsort,2));
spcsorttemp = [spcsort(:,end-plfdlg:end,:) spcsort spcsort(:,1:plfdlg,:)];
firstround = zeros(size(spcsort,1),1);
for i = 1:size(spcsort,1) %loop over cells
    
    a = find(abs(Alpha-angle(spcpol{i,1})) == min(abs(Alpha-angle(spcpol{i,1}))));
    b = permute(sum(spcsorttemp(i,a:a+plfdlg,:),2),[3 2 1]);
    if find(b,1)
        firstround(i) = find(b,1);  
    end
end
end

function [plcfield,centerfield,firstround] = dombeck(spcsort,plfdlg)
    %%    
    % lengthen spcsort to avoid errors at border
    spcsorttemp = [spcsort(:,end-plfdlg:end,:) spcsort spcsort(:,1:plfdlg,:)];
    plcfield = zeros(size(spcsort,1),size(spcsort,2),5);
    centerfield = zeros(size(spcsort,1),1);
    
    % Argument 1: Act in at least (roundperc)% rounds
    roundperc = .2;   % percentage of rounds where cell showes place activity
    for i = 1:size(spcsort,2) %loop over bins
        % read out plc related act and smooth over the placefld length
        plcfield(:,i,1) = mean(sum(spcsorttemp(:,i:i+plfdlg*2,:),2),3);
        temp = sum(spcsorttemp(:,i:i+plfdlg*2,:),2);
        temp(temp>1) = 1;
        temp = sum(temp,3)./size(spcsort,3);
        temp(temp<roundperc) = 0;
        temp(temp>0) = 1;
        plcfield(:,i,2) = temp;
    end
    
    % Arguement 2: Act at least (multfact) times higher within plc field
    multfact = 7;
    for i = 1:size(spcsort,2) %loop over bins                
        plcfield(:,i,3) = mean(plcfield(:,[1:i-plfdlg i+plfdlg:end],1),2);
    end
    temp = multfact * plcfield(:,:,3)<plcfield(:,:,1);
    plcfield(:,:,3) = temp;
    % Combine Argument 1 & 2 for potential plc fields
    plcfield(:,:,4) = plcfield(:,:,1).*plcfield(:,:,2).*plcfield(:,:,3);
    
    % Argument 3: Place fields need to have minimal length
    for i = 1:size(spcsort,1) %loop over cells
        temp = find(plcfield(i,:,4));
        for j = 1:size(temp,2) %loop found fields
            % look for start point
            if temp(j) == 1 || temp(j) == 2
                jj = temp(1);
            end
            if temp(j)>2 && plcfield(i,temp(j)-1,4) == 0 && plcfield(i,temp(j)-2,4) == 0 %&& plcfield(i,temp(j)+1,4) ~= 1
                jj = j;
            end
            % look for the end of the field and check length
            if temp(j)<size(plcfield(i,:,4),2)-1 && plcfield(i,temp(j)+1,4) == 0 && plcfield(i,temp(j)+2,4) == 0
                if j-jj>= plfdlg-2
                    plcfield(i,temp(jj):temp(j),5) = plcfield(i,temp(jj):temp(j),4);
                end
            % Check for place fields close to the border 
            elseif temp(j) == size(plcfield(i,:,4),2)-1
                if j-jj>= plfdlg-2
                    plcfield(i,temp(jj):temp(j),5) = plcfield(i,temp(jj):temp(j),4);
                end
            % Check for place fields at the border of the rounds
            elseif temp(j) == size(plcfield(i,:,4),2) 
                if j-jj>= plfdlg-2
                    plcfield(i,temp(jj):temp(j),5) = plcfield(i,temp(jj):temp(j),4);
                end
                if temp(1) == 1
                    jjj = find(plcfield(i,:,4)==0,1)-1;
                    if size(plcfield(i,:,4),2)-jj + jjj >= plfdlg-2
                        plcfield(i,temp(jj):end,5) = plcfield(i,temp(jj):end,4);
                        plcfield(i,1:jjj,5) = plcfield(i,1:jjj,4);
                    end
                end
            end            
        end
        temp = find(plcfield(i,:,5));
        for j = 1:size(temp,2)-1 %loop found fields to fill small gaps
            if temp(j)+2==temp(j+1)
                plcfield(i,temp(j)+1,5) = plcfield(i,temp(j)+1,1);
            end
        end
        temp = find(plcfield(i,:,5));
        if ~isempty(temp)
            centerfield(i) = mean(temp);
        end
    end    
    % Read out the first round when activity happens in plcfield 
    firstround = zeros(size(spcsort,3),1);
    for i = 1:size(spcsort,1) %loop over cells
        if find(plcfield(i,:,5)>0,1)
            a = permute(sum(spcsort(i,plcfield(i,:,5)>0,:),2),[3 2 1]);
            if isempty(find(a,1)); a = permute(sum(spcsort(i,plcfield(i,:,1)>0,:),2),[3 2 1]);end
            firstround(i) = find(a,1);  
        end
    end
end

function In = skaggs(spcsort,numbin,spcbin,disttemp)
    %%
    spcsort = spcsort(:,1:numbin,:);
    p_i = zeros(1,numbin); 
    for j = 1:numbin
        p_i(j) = length(find(disttemp>(j-1)*spcbin+1 & disttemp<j*spcbin));           
    end
    p_i = p_i/sum(p_i);

    lambda_i = sum(spcsort,3);
    lambda = zeros(size(lambda_i,1),1);
    In = zeros(size(lambda_i,1),1);
    for j = 1:size(lambda_i,1)
        % overall firing rate lambda (only running periods) [events/min]
        lambda(j) = sum(lambda_i(j,:).* p_i );
        In(j) = sum(lambda_i(j,lambda_i(j,:)~=0) .* log2(lambda_i(j,lambda_i(j,:)~=0)./lambda(j)) .* p_i(lambda_i(j,:)~=0) );
    end
end

function fldlength = fieldlength(scn,spcbin)
  %% Fit gauss for place field length        
        
    gaussEqn = 'a*exp(-((x-b)/(sqrt(2)*c))^2)';
    x = (1:size(scn.plcfield,2))*spcbin;
    xhlf = round(length(x)/2);
    xx = [x x+x(end)];
    
    fldlength = cell(size(scn.plcfield,1),2);
    for i = 1:size(scn.plcfield,1)
        %% smoothed actitivity in apsbins for fitting
        y = reshape(scn.plcfield(i,:,1),size(x));
        yy = [y(xhlf+1:end) y y(1:xhlf)];
        % maximum as peak strating point
        [a,b] = max(y);
        startPoints = [a (b+xhlf)*spcbin 4+spcbin];
        upperBounds = [a 3/4*xx(end) x(end)/4];
        lowerBounds = [0 1/4*xx(end) 0];
        f1 = fit(xx',yy',gaussEqn,...
            'Start', startPoints,...
            'Upper',upperBounds,...
            'Lower',lowerBounds);
        
        yy = f1(xx);
        fldlength{i,1} = yy(xhlf+1:xhlf+length(x));
        fldlength{i,2} = [round(f1.b-xhlf*spcbin) round(f1.c)];
    end
      
end

function [trace,traceErr,onsets,c,p,bin] = SpeedCorr(caim,scn)
    % speed is binned in cm/s
    speed = scn.speed*100;
    speed(speed<-20) = 0;
    speed(end+1)=speed(end);
    speed = smooth(speed,30);

    speed = speed(1:size(caim.C,2));
    speed(speed<=0) = 0;
    bin = min(speed) :max(speed)/20: max(speed);
    trace = zeros(size(caim.C,1),length(bin)-1);
    onsets = zeros(size(caim.C,1),length(bin)-1);
    traceErr = zeros(size(caim.C,1),length(bin)-1);
    for i = 1:length(bin)-1
        for j = 1:size(trace,1) 
            trace(j,i) = mean(caim.C(j,speed>bin(i)&speed<=bin(i+1))./caim.Df(j));
            traceErr(j,i) = std(caim.C(j,speed>bin(i)&speed<bin(i+1))./caim.Df(j));%/sqrt(sum(speed>bin(i)&speed<bin(i+1)));
            onsets(j,i) = sum(caim.S_bin(j,speed>bin(i)&speed<=bin(i+1)));  
        end
    end
    [c,p] = corr(bin(2:end)',trace');
    c = reshape(c,length(c),1);
    p = reshape(p,length(p),1);
end

function DataOut = ShuffleData(DataIn)
    caim = DataIn;
    S_bin = caim.S_bin;
    % shuffle everything
    for j = 1:size(S_bin,1)
        a = randperm(size(S_bin,2));
        caim.S_bin(j,:) = S_bin(j,a);
%         caim.C(j,:) = caim.C(j,a);
    end

%     shift traces randomly with respect to each other
    for j = 1:size(S_bin,1)
        a = randperm(size(S_bin,2),1);
%         caim.S_bin(j,:) = S_bin(j,[a:end 1:a-1]);
        caim.C(j,:)  = caim.C(j,[a:end 1:a-1]);
    end

    caim.shuffle = 1; 
    DataOut = caim;
end

function polarplots(scn)

%% Polar plots

%     spcpol = scn.spcpol;
%     j = 4;
%     jj = 1;
%     for i = 1:size(spcpol,1)
%         if j == 4
%             figure('position',[400   260   1000   700])
%             j = 1;
%             jj = 1;
%         end
%         if mod(i,3)==0
%             j = j+1;
%         end
%         d = (length(rnd)+1)*spcpol{i,1};
%         c = spcpol{i,3};
%         subplot(3,3,jj)
%         title(['Cell ' num2str(i)])
%         polarplot(theta,rho) 
%         hold on               
%         polarplot(c(:,1),c(:,2),'o')
%         polarplot(d,'*');
%         polarplot([0 real(-1j*log(d))],[0 abs(d)]);
%         axis off
%         hold off
%         jj = jj+1;
%     end
    
%% Shuffle plots
    
%     j = 4;
%     jj = 1;
%     for i = 1:size(scn.spcpol,1)
%         if j == 4
%             figure('position',[400   260   1000   700])
%             j = 1;
%             jj = 1;
%         end
%         if mod(i,3)==0
%             j = j+1;
%         end
%         subplot(3,3,jj)
%         title(['Cell ' num2str(i)])
%         a = sum(scn.shufield(i,2,:))/(size(scn.shufield(i,2,:),3));
% 
%         histogram(scn.shufpol(i,2,:))
%         hold on    
%         plot([scn.spcpol{i,2} scn.spcpol{i,2}],[0 150])
%         xlim([0 1])
%         hold off
%         title(['Prob of placefield ' num2str(round(100*a)) ' %']);
%         
% %         histogram(angle(scn.shufpol(i,1,:)))
% %         hold on   
% %         plot([angle(scn.spcpol{i,1}) angle(scn.spcpol{i,1})],[0 150])
% %         hold off        
%         jj = jj+1;      
%         
%     end
    
%% Place field analysis

%     plcfield = scn.plcfield;
%     figure('position',[34 561 1815 420],'color',[1 1 1]);
%     
%     subplot(1,5,1);imagesc(plcfield(:,:,1));
%     title('Mean averaged activity')
%     subplot(1,5,2);imagesc(plcfield(:,:,2));
%     title('Activity in 20 % of rounds')
%     subplot(1,5,3);imagesc(plcfield(:,:,3));
%     title('Contrast within activity')
%     subplot(1,5,4);imagesc(plcfield(:,:,4));
%     title('Putativ fields')
%     subplot(1,5,5);imagesc(plcfield(:,:,5));
%     title('Fields with right length')

%     exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\place criteria')
%% Place related firing

%     figure
%     imagesc(sum(spcsort,3))
%     colorbar('XColor',[1 1 1])
%     set(gca,'fontsize',20,...
%                 'YColor',[1 1 1],...
%                 'XColor',[1 1 1],...
%                 'Color',[0 0 0],...
%                 'XTickLabel',30:30:150)
%     set(gcf,'color',[0 0 0])

    % exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\place matrix')
    
%% Cell summary plot
    
    j = 1;
    spcpol = scn.spcpol;
    figure('position',[400   260   1000   700])
    for i = 1:size(scn.spcpol,1)       
        if j > 9
            j = 1;
            figure('position',[400   260   1000   700],'color',[1 1 1])
        end
        
        subplot(3,3,j)        
        plot(scn.plcfield(i,:,1),'color',[.5 .5 .5])
        hold on
        plot(scn.spcpol{i,4},'b')
        plot(scn.plcfield(i,:,5),'r','linewidth',2)
        title(['Cell ' num2str(scn.cellID(i,1)) ', Prob of placefield ' num2str(round(100*scn.cellID(i,3))) ' %'])
        
        subplot(3,3,j+1)
        d = (max(scn.rounds)+1)*spcpol{i,1};
        c = spcpol{i,3};       
        title(['Cell ' num2str(i)])
        polarplot(scn.theta,scn.rho) 
        hold on               
        polarplot(c(:,1),c(:,2),'o')
        polarplot(d,'*');
        polarplot([0 real(-1j*log(d))],[0 abs(d)]);
        axis off
        hold off
        
        
        subplot(3,3,j+2)
        title(['Cell ' num2str(i)])
        a = length(find(scn.shufpol(i,2,:)>scn.spcpol{i,2}))/(size(scn.shufpol(i,2,:),3));
        histogram(scn.shufpol(i,2,:),...
            'Normalization','Probability')
        hold on    
        plot([scn.spcpol{i,2} scn.spcpol{i,2}],[0 .2])
        xlim([0 1])
        hold off
        title(['Prob of stronger coding ' num2str(round(100*scn.cellID(i,6))) ' %']);
        
    
        j = j+3;      
%         exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\place field3')
    end

%% place vector order plot
    spcpol = scn.spcpol;
    [a,b] = sort((angle(cell2mat(spcpol(:,1)))));
    sigID = scn.cellID(b,:); 
    a = a(sigID(:,6)<=.05);
    sigID = sigID(sigID(:,6)<=.05,1);   
    PVorder = zeros(length(sigID),size(scn.plcfield,2));
    for i = 1:length(sigID)
%         a = find(scn.cellID(:,1)==sigID(i));
        PVorder(i,:) = scn.plcfield(scn.cellID(:,1)==sigID(i),:,1);
        PVorder(i,:) = PVorder(i,:)-min(PVorder(i,:));%/max(PVorder(i,:)-min(PVorder(i,:)));
    end
    c = find(abs(a)==min(abs(a)));
    figure;
    imagesc([PVorder(c:end,:) ; PVorder(1:c-1,:)])
    title('Sorted place vectors')
%     exportfig('C:\Users\Admin\Dropbox\Lab\Vortr�ge\place field')

%% place field order plot
    sigID = scn.cellID; 
    sigID = sigID(sigID(:,2)>0 & sigID(:,3)<=.02,:);
    [~,b] = sort(sigID(:,2));
    sigID = sigID(b,1);
    PForder = zeros(length(sigID),size(scn.plcfield,2));
    for i = 1:length(sigID)
%         a = find(scn.cellID(:,1)==sigID(i));
        PForder(i,:) = scn.plcfield(scn.cellID(:,1)==sigID(i),:,1);
        PForder(i,:) = PForder(i,:)-min(PForder(i,:));%/max(PVorder(i,:)-min(PVorder(i,:)));
    end
    figure;
    imagesc(PForder)
    title('Sorted place fields')
end


