function [] = subSystemMets(model)

% Lists metabolites involved in a certain sub system of model
%   by Gunvor Røkke, NTNU, 2021

subSysList = unique(model.subSystems);
disp([num2cell(transpose(1:length(subSysList))), subSysList])
indx = input('Index of sub sustem to be investigated: ');
disp(' ')
subSys = char(subSysList(indx));
rxns = smatch(model.subSystems, subSys, 'exact');
metList = cell.empty(0,1);
n = 1;
for j = 1:length(rxns)
    rxn = rxns(j);
    mets = find(model.S(:,rxn) ~= 0);
    for i = 1:length(mets)
        met = mets(i);
        metList(n,1) = model.metNames(met);
        n = n + 1;
    end
end
metList = unique(metList);
fprintf('Metabolites involved in sub system ''%s'':\n\n', subSys);
disp(metList)