function model = addGenes(model, rxID, geneList, andOrRelationships)

% Taking in list of genes associated with reaction
% Adding genes in geneList to model (if not already there)
% Finding index of genes in model.genes
% Filling rxnGeneMatrix for reaction in question
% Generating rule and grRule for reaction in question
%   
% Adding rule and grRule to model
% Scripts required to run addGenes:
%   - smatch.m

%   by Gunvor Røkke, NTNU, 2021
disp(' ')
rxID = char(rxID);
rxIndex = smatch(model.rxns, rxID, 'exact');
if length(rxIndex) > 1
    fprintf('Multiple reactions exist for %s\n\n', rxID)
    return
end
if ischar(geneList)
    geneList = cellstr(geneList);
end
geneIndex = zeros(length(geneList),1);
for i = 1:length(geneList)
    if isempty(smatch(model.genes, char(geneList(i)), 'exact')) % Adding gene to model.genes if it is not present in model already
        model.genes(length(model.genes)+1) = geneList(i);
    end
    if length(smatch(model.genes, char(geneList(i)), 'exact')) == 1
        geneIndex(i) = smatch(model.genes, char(geneList(i)), 'exact');
    else
        fprintf('\nDuplicate entries found in model.genes for gene %s\n\n', char(geneList(i)))
    end
end

if ~isempty(find(model.rxnGeneMat(rxIndex,:) ~= 0))
    model.rxnGeneMat(rxIndex,:) = zeros(1,size(model.rxnGeneMat,2)); % Reseting model.rxnGeneMat entry for reaction
end
for i = 1:length(geneIndex) % Filling model.rxnGeneMat
    model.rxnGeneMat(rxIndex,geneIndex(i)) = 1;
end

% Checking if andOrRelationships is present. If not, generating

if length(geneList) > 1
    if ~exist('andOrRelationships', 'var')
        createAndOr = input('Create automatic and (a) / or (o) relationships, or define manually (m)? ', 's');
        disp(' ')
        if strcmp(createAndOr, 'a')
            andOrRelationships = cell(1,(length(geneList) - 1));
            for i = 1:length(andOrRelationships)
                andOrRelationships(i) = cellstr('&');
            end
        elseif strcmp(createAndOr, 'o')
            andOrRelationships = cell(1,(length(geneList) - 1));
            for i = 1:length(andOrRelationships)
                andOrRelationships(i) = cellstr('|');
            end
        end
    elseif exist('andOrRelationships', 'var')
        if iscell(andOrRelationships) && ~isempty(smatch(andOrRelationships, '&')) && ~isempty(smatch(andOrRelationships, '|'))
            createAndOr = 'm';
        elseif ischar(andOrRelationships)
            createAndOr = 'm';
        elseif iscell(andOrRelationships)
            createAndOr = 'y';
        end
    end
end

% Building rule and grRule

if length(geneList) > 1
    if ~strcmp(createAndOr, 'm')
        rule = sprintf('(x(%i)', geneIndex(1));
        grRule = sprintf('(%s', char(geneList(1)));
        for i = 2:length(geneList)
            rule = sprintf('%s %s x(%i)', rule, char(andOrRelationships(i-1)), geneIndex(i));
            if strcmp(char(andOrRelationships(i-1)), '&')
                grRule = sprintf('%s AND %s', grRule, char(geneList(i)));
            elseif strcmp(char(andOrRelationships(i-1)), '|')
                grRule = sprintf('%s OR %s', grRule, char(geneList(i)));
            end
        end
        rule = sprintf('%s)', rule);
        grRule = sprintf('%s)', grRule);
    end
elseif length(geneList) == 1
    if ~exist('createAndOr', 'var')
        rule = sprintf('(x(%i))', geneIndex);
        grRule = sprintf('(%s)', char(geneList));
    end
end
if exist('createAndOr', 'var')
    if strcmp(createAndOr, 'm')
        disp('Create rule and grRule manually')
        disp('Gene list:')
        disp(geneList)
        disp('Gene index')
        disp(geneIndex)
        if exist('andOrRelationships', 'var')
            fprintf('And / or relationships:\n\n')
            disp(andOrRelationships)
        end
        rule = input('Define rule: ', 's');
        grRule = input('Define grRule: ', 's');
    end
end

model.rules(rxIndex) = cellstr(rule);
model.grRules(rxIndex) = cellstr(grRule);

fprintf('\n\n*** Genetic information added for reaction %s ***\n\n', char(model.rxns(rxIndex)))
fprintf('rxnGeneMat:\n\n')
disp(find(model.rxnGeneMat(rxIndex,:) ~= 0))
fprintf('\nGenes involved in reaction according to rxnGeneMat:\n\n')
disp(model.genes(find(model.rxnGeneMat(rxIndex,:) ~= 0)))
fprintf('\nRule:\n\n')
disp(model.rules(rxIndex))
fprintf('\ngrRule\n\n')
disp(model.grRules(rxIndex))
disp(' ')
disp('***')
disp(' ')

clear andOrRelationships
