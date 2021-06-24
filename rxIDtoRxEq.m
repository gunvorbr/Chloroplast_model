function index = rxIDtoRxEq(rxns, model)

% Taking in list of reaction IDs
% Finding index of reaction(s) in model
% Printing out reaction equation(s) for reaction(s)

% Required scripts:
%   - ReactionNames.m
%   - smatch.m

%   by Gunvor Røkke, NTNU, 2021

if length(rxns) == 1
    if ischar(rxns) % If only one reaction is given as input, and reaction is char -> converting to cell
        rxns = cellstr(rxns);
    end
end

index = zeros(length(rxns),1);
rxList = ReactionNames(model);

for j = 1:length(rxns)
    rxID = char(rxns(j));
    currentIndex = smatch(model.rxns, rxID, 'exact');
    if length(currentIndex) == 1
        index(j) = currentIndex;
    else
        fprintf('Reaction %s occurs multiple times in the model. Please correct this.\n\n', rxID)
    end
end

eqSubset = rxList(index);
disp(' ')
disp(char(eqSubset))