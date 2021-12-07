function spike = inferSpike( data, S, inferSpikeControl)

% this code extracts the activity event from the data
% input: data: matrix of temporal components (K x T matrix)
%        S: matrix of spike probability (K x T matrix)
%        inferSpikeControl.method: infer spike method, default: foopsi
%        inferSpikeControl.stdThr: std threshold for "foopsi" 
%        inferSpikeControl.stdThr2: std threshold for "derivative" 
%        inferSpikeControl.lowpassCutoff: low pass filter cut off coefficient, default: set to -1 means no filtering
%        inferSpikeControl.frameRate: frame rate of the recording
%        inferSpikeControl.dynamicThr: dynamic control of std threshold, set to 1 to enable
%        inferSpikeControl.SNR: SNR of each component
%        inferSpikeControl.SNR0: if SNR<SNR0, dynamicThr is applied
%        inferSpikeControl.multiplier: thr=stdThr+(SNR0-SNR)*multiplier, for "foopsi"
%        inferSpikeControl.multiplier2: thr2=stdThr2+(SNR0-SNR)*multiplier2, for "derivative"
% output: spike: spike binary matrix

% Author: Weijian Yang

if isfield(inferSpikeControl,'method') method=inferSpikeControl.method;
else method='foopsi'; end

if isfield(inferSpikeControl,'stdThr') stdThr=inferSpikeControl.stdThr;
else stdThr=2.5; end

if isfield(inferSpikeControl,'stdThr2') stdThr2=inferSpikeControl.stdThr2;
else stdThr2=2.5; end

if isfield(inferSpikeControl,'lowpassCutoff') lowpassCutoff=inferSpikeControl.lowpassCutoff;
else lowpassCutoff=-1; end

if isfield(inferSpikeControl,'frameRate') frameRate=inferSpikeControl.frameRate;
else if strcmp(method,'derivative')||strcmp(method,'foopsi_derivative') disp('Please set up frame rate...'); return; end 
end

if isfield(inferSpikeControl,'dynamicThr') && inferSpikeControl.dynamicThr
    deltaSNR=inferSpikeControl.SNR0-inferSpikeControl.SNR;
    thr=stdThr+deltaSNR*inferSpikeControl.multiplier;
    thr2=stdThr2+deltaSNR*inferSpikeControl.multiplier2;
    index=find(thr<stdThr); thr(index)=stdThr;  stdThr=thr;
    index=find(thr2<stdThr2); thr2(index)=stdThr2;  stdThr2=thr2;
else
    stdThr=repmat(stdThr,1,size(data,1));
    stdThr2=repmat(stdThr2,1,size(data,1));
end

K=size(data,1);
spike=zeros(size(data));

if strcmp(method,'foopsi')
    for idx=1:K
        index=find(S(idx,:)>mean(S(idx,:))+std(S(idx,:))*stdThr(idx));
        spike(idx,index)=1;
    end
else if strcmp(method,'derivative') || strcmp(method,'foopsi_derivative')  
    % setup low pass filter
        options.nlfilt.value=40;		    % filter length
        options.lpass.value=lowpassCutoff;  % [Hz] low pass cut off frequency
        options.timeRes.value=1/frameRate;	% [s] sampling time
    % low pass filter
        if lowpassCutoff<0
            dataFiltered=data;
        else
            dataFiltered=lowpass_filter(data, options);
        end
    % run derivative        
        dataDif=gradient(dataFiltered,1,2);
        index=find(dataDif<0); 
        dataDif(index)=0;
    % thresholding
        for idx=1:K
            dataDif(idx,:)=dataDif(idx,:)/max(dataDif(idx,:));      % normalization
            if strcmp(method,'derivative')
                index=find(dataDif(idx,:)>(mean(dataDif(idx,:))+stdThr2(idx)*std(dataDif(idx,:))));   
                spike(idx,index)=1;
            else
                index2=find(dataDif(idx,:)>(mean(dataDif(idx,:))+stdThr2(idx)*std(dataDif(idx,:))));   
                index1=find(S(idx,:)>mean(S(idx,:))+std(S(idx,:))*stdThr(idx));
                index=intersect(index1,index2);
                spike(idx,index)=1;
            end
        end
    end
end

