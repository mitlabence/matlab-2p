function stim = RespClust(stim,response,imID,mouseIDtempgc,mouseIDtempbulk)

if nargin == 0
    stim.mouseID = [];
    stim.mouseIDbulk = [];
    stim.nums = [];
    stim.times   = []; 
    stim.speed   = [];
    stim.pupil   = [];  
    stim.sumnet  = [];
    stim.sumresp = [];    
    stim.nresp = [];
    stim.plcresp = [];
    stim.spcresp = [];
    stim.bulk = [];
    stim.stimplc = [];
    stim.effects = cell(0);
    stim.pretime = cell(0);
    stim.waittime = cell(0);
    stim.runtime = cell(0);
    stim.stimplc = [];
    return
end

if isempty(response.expID);return;end
expID = zeros(length(response.expID),1);
for i = 1:length(response.expID)
    expID(i) = find(response.expID(i)== cell2mat(mouseIDtempgc(:,2)),1);
end

mouseIDtempgc= mouseIDtempgc(expID,:);

stim.mouseID = [stim.mouseID; mouseIDtempgc];
stim.times   = [stim.times; response.times]; 
stim.speed   = [stim.speed; response.speed];
stim.pupil   = [stim.pupil; response.pupil];
stim.nums    = [stim.nums; response.nums];
stim.effects = [stim.effects; response.effects];
stim.pretime = [stim.pretime; response.pretime];
stim.waittime = [stim.waittime; response.waittime];
stim.runtime = [stim.runtime; response.runtime];
stim.stimplc = [stim.stimplc; response.stimplc];

if imID(1)    
    stim.sumnet  = [stim.sumnet; response.sumnet];
    stim.sumresp = [stim.sumresp; response.sumresp];
    stim.nresp = [stim.nresp; response.nresp];
    stim.plcresp = [stim.plcresp; response.plcresp];
    stim.spcresp = [stim.spcresp; response.spcresp];
end
if imID(2)
    expID = zeros(length(response.expID),1);
    for i = 1:length(response.expID)
        expID(i) = find(response.expID(i)== cell2mat(mouseIDtempbulk(:,2)),1);
    end
    stim.bulk = [stim.bulk; response.bulk];    
    mouseIDtempbulk = mouseIDtempbulk(expID,:);
    stim.mouseIDbulk = [stim.mouseIDbulk; mouseIDtempbulk];
end

end