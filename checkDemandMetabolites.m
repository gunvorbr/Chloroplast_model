function notProd = checkDemandMetabolites(model, metList)

% Checks that a list of metabolites can be produced by adding demand
% reactions consuming them.
% Returning list of metabolites that are not produced.
%   by Gunvor Røkke, NTNU, 2021

if iscell(metList)
    indexList = zeros(length(metList),1);
    for i = 1:length(metList)
        index = smatch(model.metNames, char(metList(i)), 'exact');
        if length(index) == 1
            indexList(i) = index;
        else
            fprintf('Warning: Something wrong with metabolite %s. Metabolite not found in model, or found more than once\n\n', char(metList(i)))
        end
    end
end
metList = indexList;
clear indexList index i

disp(' ')
notProd = double.empty(0,1);
n = 1;

fluxDisp = input('Display fluxes? (yes = 1, no = 0): ');
disp(' ')

for i = 1:length(metList)
    fprintf('*** Metabolite %s ***\n\n', char(metList(i)))
    if find(model.c ~= 0)
        warning('Reaction is already being optimized')
        return
    end
    model = demandRx(model, char(model.mets(metList(i))));
    model.c(end) = 1;
    FBA = optimizeCbModel(model);
    if FBA.f < 0.001
        fprintf('\nMetabolite %s cannot be produced by the model\n\n', char(model.mets(metList(i))))
        notProd(n,1) = metList(i);
        n = n + 1;
    elseif FBA.f > 0.001
        fprintf('\nMetabolite %s can be produced by the model :-D\n\n', char(model.mets(metList(i))))
    end
    % Displaying fluxes (if user wants to)
    if fluxDisp == 1
        FBAresult = optimizeCbModel(model);
        shortSubSys = cell(length(model.subSystems),1);
        for j = 1:length(model.subSystems)
            currentSubSys = char(model.subSystems(j));
            if length(currentSubSys) > 30
                shortSubSys(j) = cellstr(currentSubSys(1:30));
            elseif length(currentSubSys) < 30
                shortSubSys(j) = cellstr(currentSubSys);
            end
        end
        disp([model.rxns, num2cell(FBAresult.x), shortSubSys])
        clear FBAresult
        goFurther = input('Continue?');
    end
    % Deleting demand reaction
    model = removeRxns(model, model.rxns(end));
    if find(model.c ~= 0)
        warning('Something went wrong')
        return
    end
end
