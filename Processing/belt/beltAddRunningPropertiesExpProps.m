function [belt_struct, params] = beltAddRunningPropertiesExpProps(belt_struct, params)
%BELTADDRUNNINGFIELDEXPPROPS The same as beltAddRunningProperties, with the
% addition of a parameter structure params that contains the important
% parameters
% This function appends the fields "running" and 
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
running(smooth(abs(belt_struct.speed))>thres) = 1;
stepup = find(diff(running)==1);
stepdown = find(diff(running)==-1); 
for i = 1 : length(stepdown)-1
    if length(stepup) < i+1
        
    elseif stepup(i+1)-stepdown(i)<wind
        running(stepdown(i):stepup(i+1)) = 1;
    end
end
runtime = belt_struct.time(running == 1);

if isprop(params, "belt_input_thres")
    disp( "beltMatchToNikonStampsExpProps: belt_input_thres is overwritten!");
end
if isprop(params, "belt_interrunning_window")
    disp( "beltMatchToNikonStampsExpProps: belt_interrunning_window is overwritten!");
end
params.belt_input_thres = thres;
params.belt_interrunning_window = wind;

belt_struct.running = running;
belt_struct.runtime = runtime;

end

