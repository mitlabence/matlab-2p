function belt_struct = beltAddRunningProperties(belt_struct)
%BELTADDRUNNINGFIELD This function appends the fields "running" and 
%"runtime" to the struct belt_struct. It is 1 at each time point when the 
%mouse is running, and 0 if not. Assumes 100 Hz belt data.
% Input:
%   belt_struct: the belt data imported as a struct (i.e. the data
%   columns as fields). See readcaim.m from Martin, where Belt structure is
%   created from the belt matrix.
% Output:
%   belt_struct: the same belt_struct extended by "running" (0 when not
%   running, 1 when running) and "runtime" (see implementation).
%
% Code: from Martin's BeltToSCN ("define periods of running" part)

%TODO: add manual threshold, window

% Threshold for input signal
thres = 40;
% Width of interrunning period to be considered still running (in bins @ 100 Hz)
wind = 250;
running = zeros(length(belt_struct.speed),1);
running(smooth(belt_struct.speed)>thres) = 1;
stepup = find(diff(running)==1);
stepdown = find(diff(running)==-1); 
for i = 1 : length(stepdown)-1
    if length(stepup) < i+1
        
    elseif stepup(i+1)-stepdown(i)<wind
        running(stepdown(i):stepup(i+1)) = 1;
    end
end
runtime = belt_struct.time(running == 1);

belt_struct.running = running;
belt_struct.runtime = runtime;

end

