function [] = moveReactions(fromModel, toModelName, moveList)

% Generates syntax to move reactions from one model to another
% After script has run, paste syntax generated by the script on prompt to
% move reactions
% Inputs:
%   fromModel:      Model structure reactions should be moved FROM
%   toModelName:    Name of model reactions should be moved TO
%   moveList:       List of indexes of reactions to be moved

%   by Gunvor R�kke, NTNU, 2021

disp(' ')
for i = 1:length(moveList)
    
    rxID = char(fromModel.rxns(moveList(i))); % Reaksjons-ID fra vektor rxns
    
    sIndex = find(fromModel.S(:,moveList(i)));
    
    Sstring = sprintf('[%f', nonzeros(fromModel.S(sIndex(1), moveList(i)))); % Lokal st�kiometrisk matrise
    metString = sprintf('{''%s''', char(fromModel.mets(sIndex(1)))); % Celle med metabolitt-IDer
    metNameS = sprintf('{''%s''', char(fromModel.metNames(sIndex(1)))); % Celle med metabolittnavn
    KEGGstring = sprintf('{''%s''', char(fromModel.metKEGGID(sIndex(1)))); % Celle med KEGG-IDer
    
    if length(sIndex) > 1
        for j = 2:length(sIndex)
            Sstring = sprintf('%s, %f', Sstring, nonzeros(fromModel.S(sIndex(j), moveList(i))));
            metString = sprintf('%s, ''%s''', metString, char(fromModel.mets(sIndex(j))));
            metNameS = sprintf('%s, ''%s''', metNameS, char(fromModel.metNames(sIndex(j))));
            KEGGstring = sprintf('%s, ''%s''', KEGGstring, char(fromModel.metKEGGID(sIndex(j))));
        end
    end
    
    Sstring = sprintf('%s]', Sstring);
    metString = sprintf('%s}', metString);
    metNameS = sprintf('%s}', metNameS);
    KEGGstring = sprintf('%s}', KEGGstring);
    
    rev = fromModel.rev(moveList(i));
    
    lb = fromModel.lb(moveList(i));
    
    ub = fromModel.ub(moveList(i));
    
    subSystem = char(fromModel.subSystems(moveList(i)));
    
    rxName = char(fromModel.rxnNames(moveList(i)));
    
    ECnum = char(fromModel.rxnECNumbers(moveList(i)));
    
    rxRef = char(fromModel.rxnReferences(moveList(i)));
    
    geneHit = transpose(nonzeros(fromModel.rxnGeneMat(moveList(i),:)));
    if ~isempty(geneHit)
        geneCell = fromModel.genes(geneHit);
        if size(geneCell,1) == 1
            geneList = sprintf('{''%s''}', char(geneCell));
        else
            geneList = sprintf('{''%s''', char(geneCell(1)));
            for j = 2:size(geneCell,1)
                geneList = sprintf('%s, ''%s''', geneList, char(geneCell(j)));
            end
            geneList = sprintf('%s}', geneList);
        end 
    end
    
    fprintf('Reaction %i out of %i\n', i, length(moveList))
    disp(' ')
    
    fprintf('%s = addReaction(%s, ''%s'', %s, %s, %i, %i, %i, 0, ''%s'', '''', {}, %s, true);\n', toModelName, toModelName, rxID, metString, Sstring, rev, lb, ub, subSystem, metNameS)
    if ~isempty(ECnum) && ~isempty(geneHit)
        fprintf('%s = addRxInfo(%s, %s, %s, %s, ''%s'', ''%s'', ''%s'', %s);\n', toModelName, toModelName, metString, KEGGstring, metNameS, rxName, rxRef, ECnum, geneList)
    elseif ~isempty(ECnum) && isempty(geneHit)
        fprintf('%s = addRxInfo(%s, %s, %s, %s, ''%s'', ''%s'', ''%s'', {});\n', toModelName, toModelName, metString, KEGGstring, metNameS, rxName, rxRef, ECnum)
    elseif isempty(ECnum) && isempty(geneHit)
        fprintf('%s = addRxInfo(%s, %s, %s, %s, ''%s'', ''%s'', '''', {});\n', toModelName, toModelName, metString, KEGGstring, metNameS, rxName, rxRef)
    elseif isempty(ECnum) && ~isempty(geneHit)
        fprintf('%s = addRxInfo(%s, %s, %s, %s, ''%s'', ''%s'', '''', %s);\n', toModelName, toModelName, metString, KEGGstring, metNameS, rxName, rxRef, geneList)
    end
    disp(' ')
    disp(' ')
    clear ECnum rxRef geneList geneHit
end