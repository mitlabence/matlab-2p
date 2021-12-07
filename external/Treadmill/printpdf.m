function printpdf(filename,horizontal)
    %%
    
    if nargin <2
        horizontal = 0;
    end
%     if isempty(filename)
%         [filename, pathname] = uiputfile('*.pdf','Choose figure name','C:\Users\Admin\Desktop\');
%     end
    
    if ~isempty(filename)
        h = gcf;
        set(h, 'PaperUnits','centimeters');
        set(h, 'Units','centimeters');
        pos=get(h,'Position');
        if horizontal == 0
            set(h, 'PaperSize', .5*[42 59.4]);
            set(h, 'PaperSize', [pos(3) pos(4)]);
            set(h, 'PaperPositionMode', 'manual');
            set(h, 'PaperPosition',.5*[0 0 42 59.4]);
        else
            set(h, 'PaperSize', [59.4 42]);
            set(h, 'PaperSize', [pos(3) pos(4)]);
            set(h, 'PaperPositionMode', 'manual');
            set(h, 'PaperPosition',[0 0 42 59.4]);
        end
        set(h, 'PaperPosition',[0 0 pos(3) pos(4)])
        print(gcf, '-dpdf', [filename]);  
    end
    
end