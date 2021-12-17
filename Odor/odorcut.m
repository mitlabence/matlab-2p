function [Belt,caim] = odorcut(Belt,caim)
%ODORCUT Taken from Martin Pofahl's readcaim.m file to separate funcions
%into different files for better overview.

if max(Belt.tsscn)/1000/60 < 25 || ~isfield(Belt,'odor1') || isempty(find(Belt.odor1,1))
    disp('Step odorcut is skipped.');
    return
end
%%
starttime = 10; %starttime in minutes
starttime = starttime*60*1000;

startbelt = find(min(abs(Belt.time-starttime))==abs(Belt.time-starttime));
startscn = find(min(abs(Belt.tsscn-starttime))==abs(Belt.tsscn-starttime));

caim.Y = caim.Y(:,startscn:end);
caim.C = caim.C(:,startscn:end);
caim.S = caim.S(:,startscn:end);
caim.f = caim.f(:,startscn:end);
caim.thresh = caim.thresh(:,startscn:end);
caim.S_norm = caim.S_norm(:,startscn:end);
caim.S_bin = caim.S_bin(:,startscn:end);

Belt.tsscn = Belt.tsscn(startscn:end);
    Belt.tsscn = Belt.tsscn-Belt.tsscn(1);
Belt.round = Belt.round(startbelt:end);
    Belt.round = Belt.round - Belt.round(1);
Belt.speed = Belt.speed(startbelt:end);
Belt.distance = Belt.distance(startbelt:end);
    Belt.distance = Belt.distance - Belt.distance(1);
Belt.distancePR = Belt.distancePR(startbelt:end);
Belt.reflect = Belt.reflect(startbelt:end);
Belt.licking = Belt.licking(startbelt:end);
Belt.stripes = Belt.stripes(startbelt:end);
    Belt.stripes = Belt.stripes - Belt.stripes(1);
Belt.stripesPR = Belt.stripesPR(startbelt:end);
Belt.time = Belt.time(startbelt:end);
    Belt.time = Belt.time - Belt.time(1);
Belt.timePR = Belt.timePR(startbelt:end);
Belt.reward = Belt.reward(startbelt:end);
Belt.airpuff = Belt.airpuff(startbelt:end);
Belt.soundl = Belt.soundl(startbelt:end);
Belt.soundr = Belt.soundr(startbelt:end);


Belt.odor1 = Belt.odor1(startbelt:end);
Belt.odor2 = Belt.odor2(startbelt:end);
Belt.odor3 = Belt.odor3(startbelt:end);
Belt.odor4 = Belt.odor4(startbelt:end);
Belt.odor5 = Belt.odor5(startbelt:end);


if isfield(Belt,'pupil')
    Belt.pupil = Belt.pupil(startbelt:end);
end

disp(['Data has been shortened from starting minute ' num2str(starttime/1000/60) ' to a length of ' num2str(Belt.tsscn(end)/1000/60) ' min'])


end

