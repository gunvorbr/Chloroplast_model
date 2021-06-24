function notOK = rulesChecker(model)

% Assuming that grRules are correct
% Checking that rxnGeneMat and rules match grRules
%   by Gunvor Røkke, NTNU, 2021

rxOK = zeros(length(model.grRules),1);

for j = 1:length(model.grRules) % Looping through reactions
    if ~isempty(model.grRules{j}) % If gene info is present for reaction
        % Checking rxnGeneMat vs. grRule
        geneIndx = find(model.rxnGeneMat(j,:) ~= 0);
        geneList = model.genes(geneIndx); % List of genes involved in reaction
        grRule = char(model.grRules(j)); % Getting grRule
        checkVec = zeros(length(geneList),1);
        for i = 1:length(geneList) % Checking that every gene in geneList is present in grRule
            if isempty(strfind(grRule, char(geneList(i))))
                fprintf('\nGene %s should be present in grRule %i, but is not found\n\n', char(geneList(i)), j)
            else
                checkVec(i) = 1;
            end
        end
        if (length(strfind(lower(grRule), ' or ')) + length(strfind(lower(grRule), ' and '))) ~= length(geneList) - 1 % Checking that geneList contains the right amount of genes
            fprintf('\nNumber of and/or in grRule %i does not match the number of genes in geneList\n\n', j)
            disp(geneList)
            disp(grRule)
            disp(' ')
            lengthOK = 0;
        else
            lengthOK = 1;
        end
        % Checking rule
        rule = char(model.rules(j)); % getting rule
        ruleVec = zeros(length(geneIndx),1);
        for i = 1:length(geneIndx) % Checking that nonzero entries in rxnGeneMat matches rule
            if isempty(strfind(rule, num2str(geneIndx(i))))
                fprintf('\nGene nr. %i should be present in rule %i, but is not found\n\n', num2str(geneIndx(i)), j)
            else
                ruleVec(i) = 1;
            end
        end
        if (length(strfind(rule, '|')) + length(strfind(rule, '&'))) ~= length(geneIndx) - 1 % Checking that the amount of genes in rule and the amount of nonzero entries in rxnGeneMat is the same
            fprintf('\nNumber of |/& in rule %i does not match number of genes in rxnGeneMat\n\n', j)
            disp(geneIndx)
            disp(rule)
            disp(' ')
            ruleOK = 0;
        else
            ruleOK = 1;
        end
        if (lengthOK == 1) && isequal(ones(length(checkVec),1), checkVec) && (ruleOK == 1) && isequal(ones(length(ruleVec),1), ruleVec)
            rxOK(j) = 1;
        end
    else
        if isempty(model.rules{j}) && isequal(model.rxnGeneMat(j,:), zeros(1,size(model.rxnGeneMat,2)))
            rxOK(j) = 1;
        end
    end
end
notOK = find(rxOK == 0);
if ~isempty(notOK)
    disp('Reactions that are not correct:')
    disp(notOK)
    disp(' ')
end