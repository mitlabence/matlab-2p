function caim = convertComp(filename)

disp('Function convertComp was used (external code from Martin Pofahl). It is outdated, and calling is an issue. Consider debugging this!');
[filename1, pathname] = uigetfile('*.mat','Select Component file','Z:\imaging Negar & Martin\M103\');
if ~isempty(filename1)
    load([pathname filename1]);

    caim.A = A_act;
    caim.b = b;
    caim.C = C_act;
    caim.cID = cID;
    caim.Cn = Cn;
    caim.Df = Df_act;
    caim.f = f;
    caim.options = options;
    caim.S = S_act;
    caim.Y = Y_act;

    % a = find(filename == '.');
    % filename = [filename(1:a(2)+4) 'Ca.mat'];
    save([pathname filename 'Ca.mat'],...
                'caim');
else
    caim = [];
end

end
