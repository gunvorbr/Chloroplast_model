function list = findCompartmentRxns(model)

% Lists reactions in certain compartment(s). Includes exchange reactions.
% Output: list (5 columns)
%   Col 1: rxn
%   Col 2: rxnName
%   Col 3: subSystem
%   Col 4: rxnList
%   Col 5: rxList

%   by Gunvor Røkke, NTNU, 2021

compartments = cell(1);
fillInCompartment = 1;
c = 1;

while fillInCompartment == 1
    compartments(c,1) = cellstr(input('Compartment symbol? ', 's'));
    c = c + 1;
    disp(' ')
    fillInCompartment = input('Add more compartments? (Yes = 1, No = 0): ');
    disp(' ')
    disp(' ')
end
clear c fillInCompartment

rxnList = printRxnFormula(model);
rxList = ReactionNames(model);

list = cell(1,5);
n = 1;

for i = 1:length(rxnList)
    currentStr = char(rxnList(i));
    if length(compartments) == 1
        if ~isempty(strfind(currentStr, sprintf('[%s]', char(compartments(1)))))
            list(n,1) = model.rxns(i);
            list(n,2) = model.rxnNames(i);
            list(n,3) = model.subSystems(i);
            list(n,4) = rxnList(i);
            list(n,5) = rxList(i);
            n = n + 1;
        end
    elseif length(compartments) == 2
        if ~isempty(strfind(currentStr, sprintf('[%s]', char(compartments(1))))) || ~isempty(strfind(currentStr, sprintf('[%s]', char(compartments(2)))))
            list(n,1) = model.rxns(i);
            list(n,2) = model.rxnNames(i);
            list(n,3) = model.subSystems(i);
            list(n,4) = rxnList(i);
            list(n,5) = rxList(i);
            n = n + 1;
        end
    elseif length(comparments) == 3
        if ~isempty(strfind(currentStr, sprintf('[%s]', char(compartments(1))))) || ~isempty(strfind(currentStr, sprintf('[%s]', char(compartments(2))))) || ~isempty(strfind(currentStr, sprintf('[%s]', char(compartments(3)))))
            list(n,1) = model.rxns(i);
            list(n,2) = model.rxnNames(i);
            list(n,3) = model.subSystems(i);
            list(n,4) = rxnList(i);
            list(n,5) = rxList(i);
            n = n + 1;
        end
    end
end