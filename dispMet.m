function [] = dispMet(model, metNr)

% Displays info about metabolite in a model
% Inputs: model, metNr (metabolite number in array model.mets)
%   by Gunvor Røkke, NTNU, 2021

disp(' ')
disp(model.mets(metNr))
disp(model.metNames(metNr))
disp(model.metKEGGID(metNr))