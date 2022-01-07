function belt_struct = beltMatrixToStruct(belt)
%BELTMATRIXTOSTRUCT Convert belt raw (matrix) data to struct (with fields).
% Code from Martin's readcaim.m function.

%% create belt output in struct variable for convinient handling

belt_struct = struct;
%belt_struct.tsscn = tsscn; %removed to be consistent with function naming
belt_struct.round = belt(:,1) - belt(1,1);
belt_struct.speed = belt(:,2);
belt_struct.distance = belt(:,3) - belt(1,3);
belt_struct.distancePR = belt(:,4);
belt_struct.reflect = belt(:,5);
belt_struct.licking = belt(:,6);
belt_struct.stripes = belt(:,7) - belt(1,7);
belt_struct.stripesPR = belt(:,8);
belt_struct.time = belt(:,9);
belt_struct.timePR = belt(:,10);
belt_struct.reward = belt(:,11);
belt_struct.airpuff = belt(:,12);
belt_struct.soundl = belt(:,13);
belt_struct.soundr = belt(:,14);

if size(belt,2)>14
    belt_struct.odor1 = belt(:,15);
    belt_struct.odor2 = belt(:,16);
    belt_struct.odor3 = belt(:,17);
    belt_struct.odor4 = belt(:,18);
    belt_struct.odor5 = belt(:,19);
end

if size(belt,2)>19
    belt_struct.pupil = belt(:,20);
end

end

