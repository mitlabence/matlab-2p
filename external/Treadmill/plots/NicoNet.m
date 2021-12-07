
%%
Scont = [];
Ska = [];
num = 7;
% x =scn.tsscn/1000/60;
figure('renderer','painters')
for i = 1:length(num)
    
    fullCAIM = 1:size(CAIM,2);
    emptyCAIM = [];
    for j = 1:length(fullCAIM)
        if isempty(CAIM(num(i),fullCAIM(j)).S) %%|| isempty(CAIM(i,fullCAIM(j)).PCAout) || sum(CAIM(i,k).network.netev)<10 || CAIM(i,k).behave.numrounds<5
            emptyCAIM = [emptyCAIM fullCAIM(j)];
        end
    end
    
     
    hold on
    fullCAIM(emptyCAIM) = [];
   
    for j = 1:length(fullCAIM)    
        k = fullCAIM(j); 
        S = sum(CAIM(num(i),k).S,1);
        if k <6
            plot(cumsum(histcounts(S,0:20,'normalization','probability')),'b')
            Scont = [Scont S];
        else
            plot(cumsum(histcounts(S,0:20,'normalization','probability')),'r')
            Ska = [Ska S];
        end
    end
    xlabel('Synchronous events')
    ylabel('cum prob')
end

Scont = Scont(Scont>0);
Ska = Ska(Ska>0);

%%
figure('renderer','painters')
hold on
plot(cumsum(histcounts(Scont,0:20,'normalization','probability')))
plot(cumsum(histcounts(Ska,0:20,'normalization','probability')))
xlabel('Synchronous events')
ylabel('cum prob')
legend({'Control' 'KA'})
legend('boxoff')

% plot(cumsum(histcounts(Scont,0:.01:1,'normalization','probability')))
% plot(cumsum(histcounts(Ska,0:.01:1,'normalization','probability')))