function [nikon_time_stamps] = nikonStampsCorrectArtifacts(nikon_time_stamps)
%NIKONSTAMPSCORRECTARTIFACTS Given the nikon time stamps data, this
%function corrects some recording artifacts.
% nikon_time_stamps: the result of 
%       stamps = importdata(.../nikon_time_stamps_nik.txt);
%       stamps = stamps.data;
% the code was taken from Martin's function: readcaim.m lines 56-69

% Delete Artifact that sometimes occurs...
if isnan(nikon_time_stamps(end,end))  
    nikon_time_stamps = nikon_time_stamps(1:end-1,:);
end

% if realtime correction ist currupted, figure it out here
if nikon_time_stamps(2,2) == 0.1
    disp('Realtime correction is NIS Elements has been corrupted')
    nikon_time_stamps = nikon_time_stamps(:,3)*1000;
else
    nikon_time_stamps = nikon_time_stamps(:,2)*1000;
end

disp(['NIS Elements recorded frames: ' num2str(length(nikon_time_stamps))]);

end

