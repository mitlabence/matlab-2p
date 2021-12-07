function plotresp(input)
% resp = scn.airpuff.resp;
% times = scn.airpuff.times;
% meanresp = scn.airpuff.meanresp;

resp = input.resp;
times = input.times;
meanresp = input.meanresp;

numcells = length(resp(:,1));
numstim = length(resp(1,:));

%%
close all
for i = 1:numcells
    figure
    for j = 1:numstim
        plot(times{i,j},resp{i,j})
        hold on 
    end
    plot(times{i,1},meanresp(i,:),'linewidth',2)
    hold off
end

% %%
% for i = 1:numcells
%     figure
%     plot(meanresp(i,:))
% end
end