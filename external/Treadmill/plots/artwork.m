% Artwork
% Artwork
%  load('Z:\imaging Negar & Martin\M194\Analysis\Cues\preprocessed\M194.260717.1324pro.mat')
%  load('Z:\imaging Negar & Martin\M194\Analysis\Baseline\preprocessed\M194.210717.1312pro.mat')

C = caim.C;
for i = 1:size(C,1)
    C(i,:) = C(i,:)/caim.Df(i);
%      C(i,:) = zscore(C(i,:));   
end
[a,~,~,~,explained] = pca(C);

speed = scn.speed(run)*100;
speed(speed<0) = 0;
speed(end+1)=speed(end);
speed = smooth(speed,30);
% speed = speed-min(speed);
maxspeed = max(speed);
speed = round(speed/max(speed).*1500);
%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 2.5*8.9 2.5*8.9]],...
    'PaperUnits','centimeters',...
    'PaperSize', [2.5*8.9 2.5*8.9])
hold on

%%
% run = find(scn.running==1);
run = find(scn.speed>0);
space = round(scn.distance(run));
% space = round(scn.totdist(run));
space(space==0) = 1;

mycolormap = hsv(1500);

for j = 1:10
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot([a(run(i),j) a(run(i+1),j)],[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(space(i),:))
        end    
    end
end

%%

mycolormap = jet(10);
ofst = .3;
for j = 1:10
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot([a(run(i),j) a(run(i+1),j)],ofst+[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(j,:))
        end    
    end
end


%%

mycolormap = jet(10);
ofst = .3;
for j = 1:10
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot(ofst+[a(run(i),j) a(run(i+1),j)],[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(j,:))
        end    
    end
end

%%
% run = find(scn.running==1);
run = find(scn.speed>0);
% space = round(scn.distance(run));
space = round(scn.totdist(run));
space(space==0) = 1;


mycolormap = jet(max(space));

for j = 1:10
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot(ofst+[a(run(i),j) a(run(i+1),j)],ofst+[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(space(i),:))
        end    
    end
end

%%
ylim = [-.15 .45];
xlim = [-.15 .45];
axis off

print(gcf, '-dpdf','PCA artwork')
close all

% Artwork
%  load('Z:\imaging Negar & Martin\M194\Analysis\Cues\preprocessed\M194.260717.1324pro.mat')
%  load('Z:\imaging Negar & Martin\M194\Analysis\Baseline\preprocessed\M194.210717.1312pro.mat')
%%

C = caim.C;
for i = 1:size(C,1)
    C(i,:) = C(i,:)/caim.Df(i);
%     C(i,:) = zscore(C(i,:));   
end
C = smoother(C, 5,1);
[a,~,~,~,explained] = pca(C);
numcomp = 7;

%%
figure('color',[1 1 1],...
    'renderer','painters',...
    'visible','on',...
    'Units','centimeters',...
    'position',[10 2 [ 60 60]],...
    'PaperUnits','centimeters',...
    'PaperSize', [60 60])
hold on

%%

run = find(scn.speed>0);
space = round(scn.distance(run));
% space = round(scn.totdist(run));
space(space==0) = 1;

mycolormap = hsv(1500);

for j = 1:numcomp
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot([a(run(i),j) a(run(i+1),j)],[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(space(i),:))
        end    
    end
end

%%

mycolormap = jet(numcomp);
ofst = .3;
for jj = 1:numcomp
    j = numcomp+1-jj;
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot([a(run(i),j) a(run(i+1),j)],ofst+[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(j,:))
        end    
    end
end

%%
speed = scn.speed(run)*100;
speed(speed<0) = 0;
speed(end+1)=speed(end);
speed = smooth(speed,30);
% speed = speed-min(speed);
maxspeed = max(speed);
speed = round(speed/max(speed).*1500);

mycolormap = jet(1500);
ofst = .3;
for j = 1:numcomp
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot(ofst+[a(run(i),j) a(run(i+1),j)],[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(speed(i),:))
        end    
    end
end

%%

% mycolormap = jet(10);
% ofst = .3;
% for j = 1:10
%     for i = 1:length(run)-1
%         if diff(run(i:i+1))<5
%             plot(ofst+[a(run(i),j) a(run(i+1),j)],[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(j,:))
%         end    
%     end
% end

%%
% run = find(scn.running==1);
run = find(scn.speed>0);
% space = round(scn.distance(run));
space = round(scn.totdist(run));
space(space==0) = 1;


mycolormap = jet(max(space));

for j = 1:numcomp
    for i = 1:length(run)-1
        if diff(run(i:i+1))<5
            plot(ofst+[a(run(i),j) a(run(i+1),j)],ofst+[a(run(i),j+1) a(run(i+1),j+1)],'color',mycolormap(space(i),:))
        end    
    end
end

%%
ylim = [-.15 .45];
xlim = [-.15 .45];
axis off
%%
print(gcf, '-dpdf','PCA artwork7')
close all