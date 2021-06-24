function rxns = findRxInModel(model)

% Lets user search for reactions involving two metabolites
% Script lets user input names of the two metabolites, and searches a model
% to see if they have any reactions in common

%   by Gunvor Røkke, NTNU, 2021

%% Searching for metabolite 1 in model
disp(' ')
met1 = input('Search term metabolite 1: ', 's');
disp(' ')
if (length(met1) == 6) && ~isnan(str2double(met1(2:end)))
    KEGG = 1;
    res1 = smatch(model.metKEGGID, met1, 'exact');
else
    res1 = findPattern(model.metNames, met1);
    KEGG = 0;
end
if KEGG == 0
    if isempty(res1{1,1})
        disp('No results found :-(')
        metNr1 = 0;
    end
elseif KEGG == 1
    if isempty(res1)
        disp('No results found :-(')
        metNr1 = 0;
    end
end
if KEGG == 0
    if ~isempty(res1{1,1})
        disp(res1)
        metNr1 = input('Choose number of metabolite 1 (if metabolite not present, press 0): ');
    end
elseif KEGG == 1
    if ~isempty(res1)
        disp([res1, model.metNames(res1)])
        metNr1 = input('Choose number of metabolite 1: ');
    end
end
if metNr1 == 0
    go = input('Search again? (yes = 1, no = 0): ');
    while go == 1
        met1 = input('New search term metabolite 1: ', 's');
        disp(' ')
        res1 = findPattern(model.metNames, met1);
        disp(res1)
        if ~isempty(res1{1,1})
            metNr1 = input('Choose number of metabolite 1 (if not present, press 0): ');
            if metNr1 ~= 0
                go = 0;
            elseif metNr1 == 0
                go = input('Search again? (yes = 1, no = 0): ');
            end
        elseif isempty(res1{1,1})
            go = input('Search again? (yes = 1, no = 0): ');
        end
    end
end
disp(' ')

%% Searching for metabolite 2 in model
met2 = input('Search term metabolite 2: ', 's');
disp(' ')
if (length(met2) == 6) && ~isnan(str2double(met2(2:end)))
    KEGG = 1;
    res2 = smatch(model.metKEGGID, met2, 'exact');
else
    res2 = findPattern(model.metNames, met2);
    KEGG = 0;
end
if KEGG == 0
    if isempty(res2{1,1})
        disp('No results found :-(')
        metNr2 = 0;
    end
elseif KEGG == 1
    if isempty(res2)
        disp('No results found :-(')
        metNr2 = 0;
    end
end
if KEGG == 0
    if ~isempty(res2{1,1})
        disp(res2)
        metNr2 = input('Choose number of metabolite 2 (if metabolite not present, press 0): ');
    end
elseif KEGG == 1
    if ~isempty(res2)
        disp([res2, model.metNames(res2)])
        metNr2 = input('Choose number of metabolite 2: ');
    end
end
if metNr2 == 0
    go = input('Search again? (yes = 1, no = 0): ');
    while go == 1
        met2 = input('New search term metabolite 2: ', 's');
        disp(' ')
        res2 = findPattern(model.metNames, met2);
        disp(res2)
        if ~isempty(res2{1,1})
            metNr2 = input('Choose number of metabolite 2 (if not present, press 0): ');
            if metNr2 ~= 0
                go = 0;
            elseif metNr2 == 0
                go = input('Search again? (yes = 1, no = 0): ');
            end
        elseif isempty(res2{1,1})
            go = input('Search again? (yes = 1, no = 0): ');
        end
    end
end
disp(' ')

%% Searching for reaction link between metabolites
rxns = zeros(1,1);
n = 1;

if (metNr1 ~= 0) && (metNr2 ~= 0)
    for j = 1:size(model.S,2)
        if (model.S(metNr1,j) ~= 0) && (model.S(metNr2,j) ~= 0)
            rxns(n,1) = j;
            n = n + 1;
        end
    end
    if rxns(1,1) ~= 0
        rxList = ReactionNames(model);
        for j = 1:length(rxns)
            fprintf('Reaction %i:\n%s\n\n', rxns(j), char(rxList(rxns(j))))
        end
    elseif rxns(1,1) == 0
        disp('No reaction link between metabolites found!')
    end
else
    if (metNr1 == 0) && (metNr2 == 0)
        disp('Neither metabolite 1 nor metabolite 2 is present in model!')
    elseif (metNr1 == 0) && (metNr2 ~= 0)
        disp('Metabolite 1 not present in model!')
    elseif (metNr1 ~= 0) && (metNr2 == 0)
        disp('Metabolite 2 not present in model!')
    end
end
disp(' ')