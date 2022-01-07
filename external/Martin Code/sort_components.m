function [cID,thresh] = sort_components(caim)
disp('Function sort_components called (external). This function is undocumented and should be updated (documentation of effects, refactoring) or removed!');
%TODO: what does this function do?
%%
if isfield(caim,'raw') %TODO: which function creates field "raw"?
    disp('sort_components: caim object has field "raw", i.e. is not a CMNF object!');
    C  = caim.raw.C;
    Df = caim.raw.Df;
    S  = caim.raw.S;
    Y  = caim.raw.Y;
else %for CMNF objects, this is called
    C = caim.C;
    S = caim.S;
    Y = caim.Y;
    Df = caim.Df;
end 
%%

gaussEqn = 'a*exp(-((x-b)/(sqrt(2)*c))^2)';
thresh = zeros(size(C));
cID = zeros(size(Df,1),1);
%TODO: change these into options!
sigmamult = 2;
bins = 100;
win = 1000;
numwin = ceil(size(C,2)/win);
counts = zeros(numwin,bins);
centers = zeros(numwin,bins);
threshtemp = zeros(numwin,1);

%%

for i = 1:size(C,1)
    %%
    
    Ytemp = Y(i,:);
    for j = 1:numwin
        %%       
        if j < numwin
            tempwin = 1+(j-1)*win:j*win;
        else
            tempwin = 1+(j-1)*win:size(C,2);
        end
        
        [counts(j,:),centers(j,:)] = hist(Ytemp(tempwin),bins,'normal');
    end
    
    %%
    
    parfor j = 1:numwin
        countmax = find(counts(j,:)==max(counts(j,:)),1);
        shift = 1;
        if countmax*2<length(centers(j,:))
            x = centers(j,1:countmax*2-shift);
            y = counts(j,1:countmax*2-shift);
        else
            x = centers(j,:);
            y = counts(j,:);
        end
        startPoints = [counts(j,countmax) centers(j,countmax) 100];
        f1 = fit(x',y',gaussEqn,'Start', startPoints);
        threshtemp(j) = f1.b+sigmamult*abs(f1.c);
    
%         figure
%         hist(Ytemp(tempwin),bins,'normal')
%         hold on
%         plot(x,f1(x),'linewidth',2,'color','r')
%         plot([thresh(i) thresh(i)],[0 max(f1(x))],'linewidth',2,'color','g')
% %         grid on
%         set(gca,'Xlim',[0 1000],...
%             'fontsize',20,...
%             'YColor',[1 1 1],...
%             'XColor',[1 1 1],...
%             'Color',[0 0 0])
%         set(gcf,'color',[0 0 0])
%         hold on
%         plot(x,y)
%         hold off

    end
    

    
    %%
    for j = 1:numwin
        if j < numwin
            tempwin = 1+(j-1)*win:j*win;
        else
            tempwin = 1+(j-1)*win:size(C,2);
        end
        thresh(i,tempwin) = threshtemp(j);
    end
    
    if find(C(i,:)>thresh(i,:),1)
        cID(i) = 1;
    else       
        cID(i) = 2;
    end
    
end

%%
a = find(cID==1);
for j = 1:length(a)
    %%
    i = a(j);
    b = zeros(1,length(S(i,:)));
    b(C(i,:)>thresh(i)) = 1;
    b(1:end-2) = b(3:end);
    if isempty(find(S(i,:).*b,1))
        cID(i) = 2;
    end
end

end

