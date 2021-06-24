function problemMets = pathwayProblemMets(model, pathwayRxList)

% Identifies metabolites in a certain pathway that cannot be generated by
% the model
%   by Gunvor R�kke, NTNU, 2021

metList = double.empty(0,1);
n = 1;

for k = 1:length(pathwayRxList)
    j = pathwayRxList(k);
    for i = 1:size(model.mets,1)
        if model.S(i,j) ~= 0
            metList(n,1) = i;
            n = n + 1;
        end
    end
end
metList = unique(metList);
clear n i j k

problemMets = double.empty(0,1);
n = 1;
for i = 1:length(metList)
    if find(model.c ~= 0)
        warning('Muffins!')
        return
    end
    model = demandRx(model, char(model.mets(metList(i))));
    model.c(end) = 1;
    FBA = optimizeCbModel(model);
    if FBA.f < 0.0001
        problemMets(n,1) = metList(i);
        n = n + 1;
    end
    model = removeRxns(model, model.rxns(end));
    if find(model.c ~= 0)
        warning('Muffins!')
        return
    end
end