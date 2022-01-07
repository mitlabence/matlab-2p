function tmscn = timestampCorrect(nikon_time_stamps,belt)
%TIMESTAMPCORRECT taken from Martin Pofahl's readcaim.m file to put
%   separate functions in separate files for better overview.

%TODO: refactor it so functionality is visible!

    if find(belt(:,20),1)
        %%        
        disp('Realtime feedback from microscope did not work. Time is alighned to last pupil detection')
        puponset = find(belt(:,20));
%         pupdiff = diff(puponset);
        pupoffset = puponset(end-1);
        tmstop = belt(pupoffset,9); %timepoint of last pupil detection in belt measurement
        tmscn = flip(tmstop - nikon_time_stamps);
        if tmscn(1)<0
            pupoffset = puponset(1);
            tmstart = belt(pupoffset,9); %timepoint of first pupil detection in belt measurement
            tmscn = tmstart + nikon_time_stamps;
        end
        
        figure
        plot(belt(:,9),belt(:,20))
        hold on
%         plot([tmscn(end) tmscn(end)],[0 max(belt(:,20))])
        plot([tmscn(1) tmscn(1)],[0 max(belt(:,20))],[tmscn(end) tmscn(end)],[0 max(belt(:,20))])
        hold off
        
    end
end

