function [] = geneToRx(model)

% Displays reactions associated with a certain gene
% Script allows user to search for gene in model
%   by Gunvor Røkke, NTNU, 2021

geneSearch = 1;
while geneSearch == 1
    gene = input('Gene name: ', 's');
    disp(' ')
    res = smatch(model.genes, gene);
    disp([transpose(num2cell(1:length(res))), model.genes(res)])
    if ~isempty(res)
        geneIndex = res(input('Correct gene (hit 0 to search again): '));
        disp(' ')
        if geneIndex > 0
            geneSearch = 0;
        end
    else
        geneSearch = input('Search again? (1 = yes, 0 = no): ');
        disp(' ')
    end
end
if exist('geneIndex', 'var')
    if geneIndex > 0
        rxns = find(model.rxnGeneMat(:,geneIndex) ~= 0);
        rxList = ReactionNames(model);
        fprintf('Gene %s is involved in the following reactions:\n\n', gene)
        disp(char(rxList(rxns)))
        disp(' ')
    end
end