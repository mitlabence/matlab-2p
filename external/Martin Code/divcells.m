function caim = divcells(caim,cID)        
%%   
if isfield(caim,'raw')
    if length(caim.raw.cID)>length(cID)
        caim.raw.cID(caim.raw.cID==1) = cID;
        cID = caim.raw.cID;
    end
    caim.A = caim.raw.A;
    caim.C  = caim.raw.C;
    caim.Df = caim.raw.Df;
    caim.S  = caim.raw.S;
    caim.Y  = caim.raw.Y;
    caim.thresh = caim.raw.thresh;
end
%%
caim.raw.A = caim.A;
caim.raw.C  = caim.C;
caim.raw.Df = caim.Df;
caim.raw.S  = caim.S;
caim.raw.Y  = caim.Y;
caim.raw.cID = cID;
caim.raw.thresh = caim.thresh;

caim.A = caim.A(:,cID==1);
caim.C  = caim.C(cID==1,:);
caim.Df = caim.Df(cID==1);
caim.S  = caim.S(cID==1,:);
caim.Y  = caim.Y(cID==1,:);
caim.cID = cID(cID==1,:);
caim.thresh = caim.thresh(cID==1,:);

caim.b = caim.b;
caim.Cn = caim.Cn;
caim.f = caim.f;
caim.options = caim.options;
    
end

