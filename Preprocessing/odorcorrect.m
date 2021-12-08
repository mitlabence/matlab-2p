function belt = odorcorrect(belt)
%ODORCORRECT Taken from Martin Pofahl's readcaim.m file, to put separate
%functions in separate files for the sake of better overview.


% if isempty(find(belt(40:end-2,15:19),1))
%     return
% end

time = belt(:,9);
a = max(time);
b = find(time == a);
%%
if length(time)>b+1 && time(b +100) == 0
    belt(1:b,15:19) = belt(end-b+1:end,15:19);
    belt = belt(1:b,:);
    disp('Belt file was corrected for weird odor software artifact')   
end
%%
% plot(belt(:,9)/1000/60)
% hold on
% plot(belt(:,15)*20)
% plot(belt(:,16)*20)
% plot(belt(:,17)*30)
% plot(belt(:,18)*40)
% plot(belt(:,19)*50)
% grid on
% hold off

end