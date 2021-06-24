function list = generalExchangeFind(model)

% Identifies EVERY exchange reaction in the model, regardless of compartment
%   by Gunvor Røkke, NTNU, 2021

list = zeros(1,1);
n = 1;

for i = 1:size(model.S,2)
    involvedMets = find(model.S(:,i) ~= 0);
    metIDs = model.mets(involvedMets);
    compSymb = cell(length(metIDs),1);
    testMet = char(metIDs(1));
    if strcmp(testMet(end), ']')
        bracketStyle = 1;
    else
        bracketStyle = 0;
    end
    for j = 1:length(metIDs)
        metID = char(metIDs(j));
        if bracketStyle == 1
            compSymb(j) = cellstr(metID(end-1));
        elseif bracketStyle == 0
            compSymb(j) = cellstr(metID(end));
        end
    end
    if length(unique(compSymb)) > 1
        list(n,1) = i;
        n = n + 1;
    end
end