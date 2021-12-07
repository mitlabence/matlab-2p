function hWaitFig=waitbarTimeRemaining2(waitbarProgress, waitBarTitle)
%function waitbarTimeRemaining(waitbarProgress)
%
% Functional purpose: To dispaly a modified waitbar with values of "ealpsed
%   time" and "remaining time"
%
% Input arguments:
%   hWaitbar- a handle to matlab waitbar (it is assumed that appropriate title is added to the waibar figure.
%
%   hTic- a handle to a tic enabled at the begining of the measured processes.
%
%   waitbarProgress- value between [0,1] describing procces done vs total process that is- current index/total indexes
%
% Output Arguments: None
%
% Issues & Comments:
%
% Author and Date: Nikolay S. 2011-04-14
% Last update:     Nikolay S. 2014-08-12
%
% Updates List 
% 2014-08-12: waitbar height increased, to make all text fit the frame.
% 2012-12-18: waitbarProgress casted to double

if nargin<1
    waitbarProgress=0;
end

if nargin<2
    waitBarTitle='waitbarTimeRemaining';
    % [StackData, iStack]=dbstack(1);
    % waitBarTitle=StackData.name;
end

persistent hWaitbar;
persistent hTic;

if isempty(hWaitbar)
    % callingFunction=evalin('caller', 'mfilename');
    hWaitbar=waitbar(0, waitBarTitle, 'Name', waitBarTitle);
    set(hWaitbar, 'CloseRequestFcn',@waitbarProgressCloseFcn);

    % Increase the waibar height so the text will fit in.
    set(hWaitbar, 'Position', get(hWaitbar, 'Position')+[0, 0, 0, 5]); 
end

if isempty(hTic)
    hTic=tic;
end

if waitbarProgress>=1
    waitbarProgressCloseFcn;
%     close(hWaitbar);
%     clear('hWaitbar', 'hTic');
%     return;
end

if isempty(hTic) || isempty(hWaitbar)
    hWaitFig=[];
    return;
end

tocVal=toc(hTic)*1e-5;
% for some reason datestr fails to work with single type
waitbarProgress=max( 1e-5, double(waitbarProgress) ); % prevent devision by 0
timePassed=tocVal;
timeRemaining=(1-waitbarProgress)/waitbarProgress*tocVal;

if timePassed<1 && timeRemaining<1
    timeLegend='H:Min:Sec.mSec';
    timeFormat='HH:MM:SS.FFF';
else
    timeLegend='Day H:Min:Sec.mSec';
    timeFormat='DD HH:MM:SS.FFF';
end

waitbar( waitbarProgress, hWaitbar,...
    sprintf('Time  passed   %s  [%s].\nTime remaining %s  [%s].',...
    datestr( timePassed, timeFormat ), timeLegend,...
    datestr( timeRemaining, timeFormat ), timeLegend) );

if ~strcmpi(waitBarTitle, 'waitbarTimeRemaining') &&...
        ~strcmpi( waitBarTitle, get(hWaitbar, 'Name') )
    set(hWaitbar, 'Name', waitBarTitle);
end

    function waitbarProgressCloseFcn(~, ~) % src, evnt
        if exist('hWaitbar', 'var')==1 & ishghandle(hWaitbar)
            delete(hWaitbar);
        end
        hWaitbar=[];
        hTic=[];
        % clear('hWaitbar', 'hTic');
        % close(hWaitbar);
        % That's funny to enable...
        % beep;  % indicate end of bar
        % load('gong.mat', 'y');
        % sound( y );
    end
hWaitFig=hWaitbar;
end

