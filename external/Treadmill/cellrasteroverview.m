addpath('C:\Users\Admin\Documents\MATLAB\export fig')
mouse = [1 5 7 8  ];
exp = [3 4 5 6 7 8 9 10 11];
mouseID = {'M103', 'M155','M158','M194','M195','M224','M226','M227','M234'}; %'M229',
for i = 1:length(mouse)
    
    figure('color',[1 1 1],...
        'position',[300 100 1.5*[ 594 420]],...
        'visible','on',...
        'PaperUnits','centimeters',...
        'Units','centimeters')
    
    cellraster(CAIMcorr,mouse(i),exp);
    title(mouseID{mouse(i)});
    printpdf(mouseID{mouse(i)})
end

files = dir('*.pdf');
append_pdfs(['Cellraster.pdf'],files.name)
delete(files.name)