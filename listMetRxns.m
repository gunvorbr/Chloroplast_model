function [] = listMetRxns(model)

% Script lets user search for metabolite in model
% Script prints reactions the metabolite in question is involved in

%   by Gunvor Røkke, NTNU, 2021

rxList = ReactionNames(model);

disp(' ')
searchMet = input('Search for metabolite ...? ', 's');
disp(' ')

searchResult = findPattern(model.metNames, searchMet);

disp([transpose(num2cell(1:size(searchResult,1))), searchResult])

metIndex = input('Index of correct metabolite? ');
metNr = cell2mat(searchResult(metIndex,1));

metRxns = transpose(find(model.S(metNr,:) ~= 0));

disp(' ')
disp(char(rxList(metRxns)))
disp(' ')
disp(char(model.subSystems(metRxns)))
disp(' ')
