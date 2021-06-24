function results = knockoutImpact(model)

% Knocks out genes to assess the effect of the knockout on cell growth

% Output:
%   - results(1) = genes
%   - results(2) = growth rate under knockout
%   - results(3) = growth rate affected? (1 / 0)

%   by Gunvor Røkke, NTNU, 2021

results = cell(length(model.genes),3);

FBA = optimizeCbModel(model);

optimalGR = FBA.obj;
clear FBA

disp(' ')
for j = 1:size(model.rxnGeneMat,2)
    fprintf('Knocking out gene %i of %i\n\n', j, length(model.genes))
    rxns = find(model.rxnGeneMat(:,j) ~= 0);
    rules = model.rules(rxns);
    for i = length(rules):-1:1
        if contains(char(rules(i)), '|')
            rxns(i) = [];
        end
    end
    if ~isempty(rxns)
        bounds = [model.lb(rxns), model.ub(rxns)];
        model.lb(rxns) = zeros(length(rxns),1);
        model.ub(rxns) = zeros(length(rxns),1);
        FBA = optimizeCbModel(model);
        results(j,1) = model.genes(j);
        results(j,2) = num2cell(FBA.obj);
        if FBA.obj >= optimalGR
            results(j,3) = cellstr('0');
        elseif FBA.obj <= optimalGR
            results(j,3) = cellstr('1');
        end
        model.lb(rxns) = bounds(:,1);
        model.ub(rxns) = bounds(:,2);
        clear FBA bounds rules rxns
    else
        FBA = optimizeCbModel(model);
        results(j,1) = model.genes(j);
        results(j,2) = num2cell(FBA.obj);
        if FBA.obj >= optimalGR
            results(j,3) = cellstr('0');
        elseif FBA.obj <= optimalGR
            results(j,3) = cellstr('1');
        end
        clear FBA rules rxns
    end
end

    