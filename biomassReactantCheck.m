function noProdList = biomassReactantCheck(model)

% Loops through every metabolite used in the biomass reaction, and checks
% if the metabolite in question can be produced by the model
%   by Gunvor Røkke, NTNU, 2021

rxList = ReactionNames(model);
bioRx = smatch(model.rxns, 'bio');

if length(bioRx) == 1
    fprintf('%s\n\n', char(rxList(bioRx)))
    correct = input('Is this the correct biomass reaction? (yes = 1, no = 0): ');
    disp(' ')
    if correct == 0
        return
    end
elseif length(bioRx) > 1
    fprintf('Possible biomass reactions found in model:\n\n')
    for j = 1:length(bioRx)
        fprintf('%i: %s\n', j, char(rxList(bioRx(j))))
    end
    disp(' ')
    correct = input('Index of correct biomass reaction (0 = biomass reaction not present in list): ');
    if correct > 0
        bioRx = bioRx(correct);
    elseif correct == 0
        return
    end
elseif isempty(bioRx)
    fprintf('No possible biomass reactions found\n\n')
    return
end

reacs = find(model.S(:,bioRx) < 0);

noProdList = cell.empty(0,1);
n = 1;

optRx = find(model.c ~= 0);
model.c(optRx) = 0;

for i = 1:length(reacs)
    model = demandRx(model, char(model.mets(reacs(i))));
    model.c(end) = 1;
    FBA = optimizeCbModel(model);
    if FBA.f < 0.000001
        noProdList(n,1) = model.mets(reacs(i));
        n = n + 1;
        fprintf('\nMetabolite %s cannot be produced\n\n', char(model.metNames(reacs(i))))
    end
    demandName = sprintf('DM_%s', char(model.mets(reacs(i))));
    model = removeRxns(model, {demandName});
end