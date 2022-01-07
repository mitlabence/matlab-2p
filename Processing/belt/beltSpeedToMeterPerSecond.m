function belt_struct = beltSpeedToMeterPerSecond(belt_struct)
%BELTSPEEDTOMETERPERSECOND this function converts the belt_struct belt data
%speed field into units of m/s.
% Input:
%       belt_struct: the belt data imported as a struct (i.e. the data
%       columns named). See readcaim.m from Martin, where Belt structure is
%       created from the belt matrix.
% Output:
%       belt_struct: the same belt_struct, with speed field entries
%       modified (length unchanged).
%
% Code from "Transfer speed readout into m/s" section of Martin's 
% BeltToSCN.m

convfact = 100; % factor multiplied by LabView
belt_struct.speed = belt_struct.speed/convfact; % 
belt_struct.speed(1) = 0;
for i = 2:length(belt_struct.speed)
    belt_struct.speed(i) = belt_struct.speed(i)/(belt_struct.time(i)-belt_struct.time(i-1)); %value mm/ms
%     belt_struct.speed(i) = (belt_struct.distance(i)-belt_struct.distance(i-1))/(belt_struct.time(i)-belt_struct.time(i-1));
end

end

