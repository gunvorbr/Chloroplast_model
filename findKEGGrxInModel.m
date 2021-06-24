function possibleRxns = findKEGGrxInModel(model, KEGGrx)

% Searching for specific KEGG reaction in model.
% Tells you if reaction is present in model or not.

% Collecting KEGG IDs for metabolites involved in reaction

%   by Gunvor Røkke, NTNU, 2021

res = urlread(strcat('http://rest.kegg.jp/get/', KEGGrx));
start = strfind(res, 'EQUATION') + 12;
slutt = start;
while ~strcmp(res(slutt), sprintf('\n'))
    slutt = slutt + 1;
end
slutt = slutt - 1;
rawKEGGIDs = res(start:slutt);
rawKEGGIDs = strsplit(rawKEGGIDs);
rows = smatch(rawKEGGIDs, 'C');
KEGGs = transpose(rawKEGGIDs(rows));
clear start slutt rawKEGGIDs rows

% Removing basic metabolites (water and protons)

if smatch(KEGGs, 'C00001', 'exact')
    KEGGs(smatch(KEGGs, 'C00001', 'exact')) = [];
end
if smatch(KEGGs, 'C00080', 'exact')
    KEGGs(smatch(KEGGs, 'C00080', 'exact')) = [];
end

% Searching for metabolites in model

hits = zeros(length(KEGGs),1);
for i = 1:length(KEGGs)
    metHits = smatch(model.metKEGGID, char(KEGGs(i)), 'exact');
    if ~isempty(metHits)
        hits(i) = 1;
    end
    for k = 1:length(metHits)
        if k == 1
            rxHits = transpose(find(model.S(metHits(k),:)));
        elseif k > 1
            rxHits = [rxHits; transpose(find(model.S(metHits(k),:)))];
        end
    end
    rxHits = sort(unique(rxHits));
    if i == 1
        possibleRxns = rxHits;
    elseif i > 1
        possibleRxns = intersect(possibleRxns, rxHits);
    end
    clear metHits rxHits k
end

% Displaying result

if isempty(possibleRxns)
    if isequal(hits, ones(length(hits),1))
        fprintf('\nAll metabolites are present, but reaction is not found in model\n\n')
    else
        fprintf('\nKEGG-reaction is not found in model. Not all metabolites are present\n\n')
    end
elseif ~isempty(possibleRxns)
    % Finding KEGG reaction equation
    start = strfind(res, 'DEFINITION') + 12;
    slutt = start;
    while ~strcmp(res(slutt), sprintf('\n'))
        slutt = slutt + 1;
    end
    slutt = slutt - 1;
    rxList = ReactionNames(model);
    if length(possibleRxns) == 1
        fprintf('\nOne hit found in model: rx %i\n\n', possibleRxns)
        fprintf('KEGG reaction:       %s\n\n', res(start:slutt))
        fprintf('Reaction in model:   %s\n\n', char(rxList(possibleRxns)))
    elseif length(possibleRxns) > 1
        fprintf('\%i hits found in model\n\n', length(possibleRxns))
        fprintf('KEGG reaction:\n%s\n\n', res(start:slutt))
        for i = 1:length(possibleRxns)
            fprintf('Rx %i in model\n', possibleRxns(i))
            fprintf('%s\n\n', char(rxList(possibleRxns(i))))
        end
    end
end