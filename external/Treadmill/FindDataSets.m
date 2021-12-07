function [files,samedate,prefiles] = FindDataSets(pathname)
%%
files =dir([pathname '*.mat']);
%%
prefiles = [];
if ispc
    prefiles1 = dir([pathname 'preprocessed\*.mat']);
    for i = 1:size(prefiles1,1)
        prefiles(i).name = [prefiles1(i).folder '\' prefiles1(i).name];
    end
else
    prefiles1 = dir([pathname 'preprocessed/*.mat']);
    for i = 1:size(prefiles1,1)
        prefiles(i).name = [prefiles1(i).folder '/' prefiles1(i).name];
    end
end
%%
% return if only one file is to be red out
if length(files) == 1
    samedate{1}(1) = 1;
    samedate{2}(1) = 1;
    return
end

% read out date from filename
expdate = datetime;%(1,length(files));
for i = 1:length(files)
    a = find(files(i).name=='.');
    expdate(i) = datetime(files(i).name(a(1)+1:a(2)-1),'InputFormat','ddMMuu');
%    expdate(i) = str2double(files(i).name(a(1)+1:a(1)+2));
end
[expdate,b] = sort(expdate);
files = files(b);
numex = days(expdate(end)-expdate(1))+1;

% find files of same date
j = 1;
k = 1;
samedate = {};
for i = 1:length(files)-1
    samedate{j,1}(k) = i;
    samedate{j,2}(1) = days(expdate(i)-expdate(end))+numex;
    if expdate(i) == expdate(i+1)       
        samedate{j,1}(k+1) = i+1;
        k = k+1;
    elseif i == length(files)-1
        samederate{j+1,1}(1) = i+1;
        samedate{j+1,2}(1) = days(expdate(i+1)-expdate(end))+numex;
    else
        k = 1;
        j = j+1;
    end
end

% add dates to order in case of discarded experiment
for i = 1:length(samedate)
    expdate(1:2);
end

end