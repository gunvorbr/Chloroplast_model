function model = demandRx(model, Dmet)

% Adds extra demand reactions to a model for metabolites in list Dmet
%   by Gunvor Røkke, NTNU, 2021

rxnName = sprintf('DM_%s', Dmet);
metaboliteList = {Dmet};
stoichCoeffList = -1;
revFlag = 0;
lowerBound = 0;
upperBound = 1000;
objCoeff = 0;
subSystem = 'Demand';
grRule = '';
geneNameList = {};
systNameList = model.metNames(strmatch(Dmet, model.mets));

model = addReaction(model, rxnName, metaboliteList, stoichCoeffList, revFlag, lowerBound, upperBound, objCoeff, subSystem, grRule, geneNameList, systNameList, false);

% Fixing subSystems problem

rxNum = smatch(model.rxns, rxnName, 'exact');
model.subSystems(rxNum) = cellstr('Demand');