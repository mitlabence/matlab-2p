function [belt_struct, params] = beltCorrectLengthExpProps(belt_struct, params, actual_length_mm)
%BELTCORRECTLENGTHEXPPROPS Same as beltCorrectLength, with an extra variable
% params struct that contains some parameters about the session.
%This function modifies the belt_struct belt data so all
%its data match the actual length values (distance between stripes, ...), 
%provided as actual_length_mm.
% Input:
%       belt_struct: the belt data imported as a struct (i.e. the data
%       columns named). See readcaim.m from Martin, where Belt structure is
%       created from the belt matrix.
%       actual_length_mm: a vector containing the length of each zone (mm).
%       If left empty, a length of 1500 mm and an equal division in three
%       parts is assumed.
% Output:
%       belt_struct: the same fields as the input, with possibly modified
%       distancePR, distance fields.
%
% Code: "Normalization of belt length" from Martin's BeltToSCN.m

stripe = find(diff(belt_struct.stripes))+1;
rnd = find(diff(belt_struct.round));
numstripes = max(belt_struct.stripesPR);
% length of the belt to be normalized to
if nargin < 3 || isempty(actual_length_mm) %TODO: check if nargin < 3 is really the correct number! 3 is used in the beltCorrectLength too...
    disp("Belt length data not supported. Assuming 1.5 m 3 zones.");
    beltlgth = 1500;
    if isprop(params, "belt_length_mm")
        disp( "beltMatchToNikonStampsExpProps: belt_length_mm is overwritten!");
    end
    actual_length_mm = beltlgth/numstripes:beltlgth/numstripes:beltlgth;
    params.belt_length_mm = actual_length_mm;
else
    disp("Belt length data was supported.");
    if isprop(params, "belt_length_mm")
        disp( "beltMatchToNikonStampsExpProps: belt_length_mm is overwritten!");
    end
    params.belt_length_mm = actual_length_mm;
end
 if ~isempty(rnd)
    
    for j = belt_struct.stripesPR(1)+1:numstripes 
        i = j-belt_struct.stripesPR(1);
        if i ==1
            win = 1:stripe(i);
            if  belt_struct.distancePR(win(end)) < actual_length_mm(j)
                belt_struct.distancePR(win) = belt_struct.distancePR(win) + actual_length_mm(j) - belt_struct.distancePR(stripe(i));
            else
                belt_struct.distancePR(win) = belt_struct.distancePR(win)*actual_length_mm(j)/max(belt_struct.distancePR(win));
                belt_struct.distance(win) = belt_struct.distancePR(win)*actual_length_mm(j)/max(belt_struct.distancePR(win));
                belt_struct.distance(win(end)+1:end) = belt_struct.distance(win(end)+1:end) - (belt_struct.distance(win(end)+1)-belt_struct.distance(win(end)));
            end
        else
            win = stripe(i-1)+1:stripe(i);
            offset = belt_struct.distancePR(win(1));                               % Offset of later steps
            belt_struct.distancePR(win) = belt_struct.distancePR(win) - offset;           % Thats is substracted
            fac = (actual_length_mm(j)-actual_length_mm(j-1))/max(belt_struct.distancePR(win));
            belt_struct.distancePR(win) = fac*belt_struct.distancePR(win) + belt_struct.distancePR(win(1)-1);
        end
        
    end

    % normalize all other rounds
    for i = 1:length(rnd)-1
        for j = 1:numstripes        
            win = stripe(i*numstripes+j-1-belt_struct.stripesPR(1))+1:stripe(i*numstripes+j-belt_struct.stripesPR(1));      % define the window in one round between two stripes          
            % Correction for distance per round
            if j == 1
                offset = belt_struct.distancePR(win(1));                           % Offset of later steps
                belt_struct.distancePR(win) = belt_struct.distancePR(win) - offset;
                fac = actual_length_mm(j)/max(belt_struct.distancePR(win));                 %  round starts at zero, therefore only correction factor needed
                belt_struct.distancePR(win) = fac*belt_struct.distancePR(win);            % That is multiplied
            else
                offset = belt_struct.distancePR(win(1));                           % Offset of later steps
                belt_struct.distancePR(win) = belt_struct.distancePR(win) - offset;       % Thats is substracted
                fac = (actual_length_mm(j)-actual_length_mm(j-1))/max(belt_struct.distancePR(win));  % And the factor calculated for this new range
                belt_struct.distancePR(win) = fac*belt_struct.distancePR(win) + belt_struct.distancePR(win(1)-1); % The last poinst of the former round is added
            end
            % Corrections for summed distance 
            offset = belt_struct.distance(win(1));                             % Offset of later steps
            belt_struct.distance(win) = belt_struct.distance(win) - offset;           % Thats is substracted
            if j == 1
                fac = (actual_length_mm(j))/max(belt_struct.distance(win));    % And the factor calculated for this new range
            else
                fac = (actual_length_mm(j)-actual_length_mm(j-1))/max(belt_struct.distance(win));    % And the factor calculated for this new range
            end
            belt_struct.distance(win) = fac*belt_struct.distance(win) + belt_struct.distance(win(1)-1); % The last poinst of the former round is added
            
        end
        
    end
%%
    % normalize last round only for the summed distance
    belt_struct.distance(rnd(end)+1:end) = belt_struct.distance(rnd(end)+1:end)-belt_struct.distance(rnd(end)+1)+belt_struct.distance(rnd(end));
    belt_struct.distancePR(rnd(end)+1:end) = belt_struct.distancePR(rnd(end)+1:end)-belt_struct.distancePR(rnd(end)+1);
 end


end

