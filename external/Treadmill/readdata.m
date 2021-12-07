function [Data,Datared] = readdata(nam,options)

sframe = options.sframe;						% user input: first frame to read (optional, default 1)
crop = options.crop;
if isfield(options,'num2read')
    num2read = options.num2read;
else
    num2read=[];					% user input: how many frames to read   (optional, default until the end)
end

if strcmp(nam(end-2:end),'mat')
    load(nam)
    if exist('cdata_1','var')
        if ~isempty(num2read)
        cdata_1 = cdata_1(:,:,sframe:num2read);
        mdata_1 = mdata_1(:,:,sframe:num2read);
        end
        Y = cdata_1;
    else
        Y = cdata_2;
    end
    Y(Y<0) = 0;
    Y = uint8(round(Y));
    c = min(mdata_1,[],3);
    if sum(c(:))~=0
        Y = Y(sum(c,2)~=0,sum(c,1)~=0,:);
    else
        Y = Y(10:end-10,10:end-10,:);
    end
    pause(.1)
    clear cdata_1 mdata_1
elseif strcmp(nam(end-2:end),'tif')
    Y = bigread2(nam,sframe,num2read);
   
    if strcmp(nam(end-4),'a') && exist([nam(1:end-5) 'b.tif'],'file')
        Y1 = bigread2([nam(1:end-5) 'b.tif'],sframe,num2read);
        Y = cat(3,Y,Y1);
        clear Y1;
    end
    if strcmp(nam(end-4),'a') && exist([nam(1:end-5) 'c.tif'],'file')
        Y1 = bigread2([nam(1:end-5) 'c.tif'],1,num2read);
        Y = cat(3,Y,Y1);
        clear Y1;
    end
    if size(Y,4) == 2
        Datared = Y(1+crop(1):end-crop(2),1+crop(3):end-crop(4),:,2);
        Y = Y(1+crop(1):end-crop(2),1+crop(3):end-crop(4),:,1);
    elseif options.numchan == 2
        Datared = Y(1+crop(1):end-crop(2),1+crop(3):end-crop(4),2:2:end);
        Y = Y(1+crop(1):end-crop(2),1+crop(3):end-crop(4),1:2:end);
    else
        Y = Y(1+crop(1):end-crop(2),1+crop(3):end-crop(4),1:end);
        Datared = [];
    end    
    Y= Y-min(min(min(Y(:,:,1:100))));
    Y(Y<0)=0;
elseif strcmp(nam(end-2:end),'nd2')
    
    [Y,Datared] = nd2read(nam,sframe,num2read);
    Datared = double(Datared);
else
    error('Unknown data type')
end


Data = double(Y);

end
                              
                                         




