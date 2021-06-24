% addRxInfo.m
% Adds KEGGids in the metKEGGID vector after adding a reaction.
% Also corrects metNames, adds rxnNames, rxnReferences, ECnum (if present),
% and genetic information.
%   by Gunvor Røkke, NTNU, 2021

function model = addRxInfo(model, metList, KEGGlist, nameList, RxName, RxRef, ECnum, geneList)

metPositions = zeros(length(metList),1);

for i = 1:length(metList) % Finding positions of metabolites in met vector
    if length(smatch(model.mets, char(metList(i)), 'exact')) == 1
        metPositions(i) = smatch(model.mets, char(metList(i)), 'exact');
    else
        warning('Muffins!')
    end
end

% Changing metNames and metKEGGIDs for metabolites

for i = 1:length(metPositions)
    model.metNames(metPositions(i)) = nameList(i);
    model.metKEGGID(metPositions(i)) = KEGGlist(i);
end

% Searching for correct reaction

rxHit = 0;
n = 1;
for j = 1:size(model.S,2) % Looping through reactions
    hit = zeros(size(model.S,1),1);
    for i = 1:size(model.S,1)
        if ~isempty(find(metPositions == i)) % Row should contain S-coefficient
            if ~isempty(nonzeros(model.S(i,j)))
                hit(i) = 1;
            end
        elseif isempty(find(metPositions == i)) % Row should be empty
            if isempty(nonzeros(model.S(i,j)))
                hit(i) = 1;
            end
        end
    end
    if isequal(hit, ones(size(model.S,1),1))
        rxHit(n,1) = j;
        n = n + 1;
    end
end
if length(rxHit) > 1
    warning('Muffins! Dublettreaksjon kan eksistere i modell. Sjekk nærmere!')
    return
elseif rxHit == 0
    warning('Muffins! Reaksjonsnummer ikke identifisert!')
else
    model.rxnNames(rxHit,1) = cellstr(RxName); % Adding reaction name
    model.rxnReferences(rxHit,1) = cellstr(RxRef); % Adding reaction reference for reaction
    model.rxnECNumbers(rxHit,1) = cellstr(ECnum); % Adding EC number for reaction
    model.subSystems(rxHit,1) = model.subSystems{rxHit}; % Fixing subSystem
end

%% Handling gene information / Generating gene syntax

if ~isempty(geneList)
    if length(geneList) > 1
        for l = 1:(length(geneList) - 1)
            % Finding index of gene, and adding gene to model.genes if not already there
            if ~isempty(smatch(model.genes, char(geneList(l)), 'exact'))
                if length(smatch(model.genes, char(geneList(l)), 'exact')) == 1
                    geneNr = smatch(model.genes, char(geneList(l)), 'exact');
                else
                    warning('Gene %s found in model.genes multiple times!', char(geneList(l)))
                    disp(' ')
                end
            else
                geneNr = length(model.genes) + 1;
                model.genes(geneNr) = cellstr(geneList(l));
            end
            % Adding gene to model.rxnGeneMat
            model.rxnGeneMat(rxHit,geneNr) = 1;
            % Building syntax for model.rules and model.grRules
            if l == 1
                rule = sprintf('(x(%i)', geneNr);
                grRule = sprintf('(%s)', char(geneList(l)));
            elseif l > 1
                if strcmp(char(geneList(end)), 'o')
                    rule = sprintf('%s | x(%i)', rule, geneNr);
                    grRule = sprintf('%s OR (%s)', grRule, char(geneList(l)));
                elseif strcmp(char(geneList(end)), 'a')
                    rule = sprintf('%s & x(%i)', rule, geneNr);
                    grRule = sprintf('%s AND (%s)', grRule, char(geneList(l)));
                end
            end
        end
        rule = sprintf('%s)', rule);
        model.rules(rxHit) = cellstr(rule);
        model.grRules(rxHit) = cellstr(grRule);
    elseif length(geneList) == 1
        if ~isempty(smatch(model.genes, char(geneList), 'exact'))
            if length(smatch(model.genes, char(geneList), 'exact')) == 1
                geneNr = smatch(model.genes, char(geneList), 'exact');
            else
                warning('Gene %s found in model.genes multiple times!', char(geneList))
                disp(' ')
            end
        else
            geneNr = length(model.genes) + 1;
            model.genes(geneNr) = geneList;
        end
        % Adding gene to model.rxnGeneMat
        model.rxnGeneMat(rxHit,geneNr) = 1;
        rule = sprintf('(x(%i))', geneNr);
        grRule = sprintf('(%s)', char(geneList));
        model.rules(rxHit) = cellstr(rule);
        model.grRules(rxHit) = cellstr(grRule);
    end 
else
    if rxHit == length(model.rxns)
        model.rxnGeneMat(rxHit,:) = zeros(1,length(model.genes));
        model.rules(rxHit) = cellstr('');
        model.grRules(rxHit) = cellstr('');
    end
end

disp(' ')
rxList = ReactionNames(model);
disp(char(rxList(rxHit)))
disp(' ')
fprintf('Rx ID:   %s\n', char(model.rxns(rxHit)))
fprintf('Rx name: %s\n', char(model.rxnNames(rxHit)))
fprintf('Rx ref:  %s\n', char(model.rxnReferences(rxHit)))
if ~isempty(model.rxnECNumbers{rxHit})
    fprintf('EC num:  %s\n', char(model.rxnECNumbers(rxHit)))
end
disp(' ')
fprintf('Genetic information for model:\n\n')
fprintf('Gene(s):\n')
if length(geneList) > 1
    for l = 1:(length(geneList) - 1)
        fprintf('%s\n', char(geneList(l)))
    end
else
    fprintf('%s\n', char(geneList))
end
fprintf('grRule: %s\n', char(model.grRules(rxHit)))
fprintf('rule:   %s\n\n', char(model.rules(rxHit)))
