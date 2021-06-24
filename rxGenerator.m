function [] = rxGenerator(model)

% Generates syntax for adding extra reactions to model manually
% Lets user search for reaction metabolites in model, and add new
% metabolites if one or more reaction metabolites are not already present
% in model

%   by Gunvor Røkke, NTNU, 2021

go = 1;

metList = cell(1,1);
metNames = cell(1,1);
KEGGlist = cell(1,1);
Slist = zeros(1,1);
m = 1;

disp(' ')
while go == 1
    quitWhile = 0;
    add = 0;
    while quitWhile == 0
        met = input('Metabolite? ', 's');
        disp(' ')
        res = findPattern(model.metNames, met); % Search for metabolite
        if isempty(res{1,1}) % Metabolite not found
            disp('Metabolite not found :-(')
            disp(' ')
            again = input('Search again? (yes = 1, no = 0): ');
            disp(' ')
            if again == 0
                add = input('Add new metabolite? (yes = 1, no = 0): ');
                quitWhile = 1;
            end
        else % Metabolite found
            disp([transpose(num2cell(1:size(res,1))), res])
            pick = input('Correct metabolite # (if not in list, choose 0): ');
            if pick ~= 0
                nr = cell2mat(res(pick,1));
                metList(m,1) = model.mets(nr);
                metNames(m,1) = model.metNames(nr);
                KEGGlist(m,1) = model.metKEGGID(nr);
                quitWhile = 1;
            elseif pick == 0
                again = input('Search again? (yes = 1, no = 0): ');
                if again == 0
                    newCompMet = input('Copy already existing metabolite to new compartment (1) or Add new metabolite (2)? ');
                    if newCompMet == 1
                        disp(' ')
                        disp([transpose(num2cell(1:size(res,1))), res])
                        pick = input('Correct metabolite # (if not in list, choose 0): ');
                        if pick ~= 0
                            nr = cell2mat(res(pick,1));
                            newComp = input('New compartment symbol: ', 's');
                            oID = char(model.mets(nr));
                            oName = char(model.metNames(nr));
                            metList(m,1) = cellstr(sprintf('%s[%s]', oID(1:end-3), newComp));
                            metNames(m,1) = cellstr(sprintf('%s[%s]', oName(1:end-3), newComp));
                            KEGGlist(m,1) = model.metKEGGID(nr);
                        end
                        quitWhile = 1;
                    elseif newCompMet == 2
                        add = 1;
                        quitWhile = 1;
                        disp(' ')
                    elseif newCompMet == 0
                        return
                    end
                end
            end
        end
    end
    if add == 1 % Add metabolite to model
        quitWhile = 0;
        while quitWhile == 0
            metList(m,1) = cellstr(input('Pick metabolite ID (remember compartment!): ', 's'));
            disp(' ')
            check = smatch(model.mets, char(metList(m,1)), 'exact');
            if ~isempty(check) % Metabolite ID is already in model. Display warning
                disp('Metabolite ID aldready exists in model. Please try again')
                disp(' ')
            else
                quitWhile = 1;
                disp('ID ok')
                disp(' ')
            end
        end
        metNames(m,1) = cellstr(input('Pick metabolite name (remember compartment): ', 's'));
        KEGGsearch = input('Search KEGG to find KEGG-ID of metabolite? (yes = 1, no = 0): ');
        disp(' ')
        if KEGGsearch == 1
            searchTerm = input('Search term: ', 's');
            disp(' ')
            KmetSearch(searchTerm)
            disp(' ')
        end
        KEGGlist(m,1) = cellstr(input('KEGG ID for metabolite: ', 's'));
        disp(' ')
        if smatch(model.metKEGGID, char(KEGGlist(m,1)), 'exact')
            warning('Metabolite seems to be exist in model already. Please confirm that this is correct')
            confirm = input('Ok? (1 = yes, 0 = no): ');
            disp(' ')
            if confirm == 0
                KEGGlist(m,1) = cellstr(input('KEGG ID for metabolite: ', 's'));
            end
        end
    end
    Slist(m,1) = input('Stoichiometric coefficient for metabolite: ');
    m = m + 1;
    go = input('Continue searching for metabolites? (1 = yes, 0 = no): ');
    disp(' ')
end

samList = [num2cell(Slist), metList, metNames, KEGGlist];
m = 1;
p = 1;

for j = 1:size(samList,1)
    if cell2mat(samList(j,1)) < 0
        mList(m,:) = samList(j,:);
        m = m + 1;
    elseif cell2mat(samList(j,1)) > 0
        pList(p,:) = samList(j,:);
        p = p + 1;
    end
end

if (exist('mList', 'var') && exist('pList', 'var')) && (size(mList,1) + size(pList,1) == size(samList,1))
    samList = [mList;pList];
elseif exist('mList', 'var') && (size(mList,1) == size(samList,1))
    samList = mList;
elseif exist('pList', 'var') && (size(pList,1) == size(samList,1))
    samList = pList;
end

Slist = cell2mat(samList(:,1));
metList = samList(:,2);
metNames = samList(:,3);
KEGGlist = samList(:,4);

if size(samList,1) > 1
    metString = sprintf('{''%s'', ', char(metList(1)));
    metNameS = sprintf('{''%s'', ', char(metNames(1)));
    KEGGstring = sprintf('{''%s'', ', char(KEGGlist(1)));
    Sstring = sprintf('[%s, ', num2str(Slist(1)));
    
    for j = 2:(length(metList)-1)
        metString = sprintf('%s''%s'', ', metString, char(metList(j)));
        metNameS = sprintf('%s''%s'', ', metNameS, char(metNames(j)));
        KEGGstring = sprintf('%s''%s'', ', KEGGstring, char(KEGGlist(j)));
        Sstring = sprintf('%s%s, ', Sstring, num2str(Slist(j)));
    end
    
    metString = sprintf('%s''%s''}', metString, char(metList(end)));
    metNameS = sprintf('%s''%s''}', metNameS, char(metNames(end)));
    KEGGstring = sprintf('%s''%s''}', KEGGstring, char(KEGGlist(end)));
    Sstring = sprintf('%s%s]', Sstring, num2str(Slist(end)));
elseif size(samList,1) == 1
    metString = sprintf('{''%s''}', char(metList(1)));
    metNameS = sprintf('{''%s''}', char(metNames(1)));
    KEGGstring = sprintf('{''%s''}', char(KEGGlist(1)));
    Sstring = sprintf('%s', num2str(Slist(1)));
end

disp(KEGGstring)
disp(' ')

search = input('Search for reaction in KEGG? (1 = yes / 0 = no): ');

if search == 0
    flag = 0;
    while flag == 0
        rxnName = input('Reacion ID? ', 's');
        if smatch(model.rxns, rxnName, 'exact')
            fprintf('Reaction ID is already found in model. Please choose another one.\n\n')
        else
            flag = 1;
        end
    end
    reactionName = input('Reaction name? ', 's');
elseif search == 1
    rxns = KEGGrxGet(model);
    disp(' ')
    disp(' ')
    rxnNum = input('Index of correct reaction? ');
    disp(' ')
    disp(' ')
    rxnName = char(rxns(rxnNum));
    url = 'http://rest.kegg.jp/get/';
    res = urlread(strcat(url, rxnName));
    start = strfind(res, 'NAME') + 12;
    slutt = strfind(res, 'DEFINITION') - 2;
    disp(res(start:slutt))
    reactionName = input('Define reaction name: ', 's');
end

revFlag = input('Reversibility? ');
disp(' ')
mode = input('Run reaction in normal mode (1) or special mode (0)? ');
disp(' ')
if (mode == 1) && (revFlag == 1)
    upperBound = num2str(1000);
    lowerBound = num2str(-1000);
elseif (mode == 1) && (revFlag == 0)
    upperBound = num2str(1000);
    lowerBound = num2str(0);
elseif mode == 0
    upperBound = num2str(0);
    lowerBound = num2str(0);
end
revFlag = num2str(revFlag);
objCoeff = num2str(0);
dispSub = input('Display sub systems for inspiration? (yes = 1, no = 0): ');
disp(' ')
if dispSub == 1
    disp(char(unique(model.subSystems)))
    disp(' ')
end
subSystem = input('Subsystem? ', 's');

if search == 1
    start = strfind(res, 'ENZYME');
    if ~isempty(start)
        ECs = res(start+12:start+50);
        ECvec = strsplit(ECs);
        ECnum = char(ECvec(1));
    elseif isempty(start)
        disp(' ')
        disp('EC number not found :-(')
        ECnum = '';
        disp(' ')
    end
elseif search == 0
    ECnum = input('EC number? ', 's');
    disp(' ')
end

rxRef = input('Reaction reference? (hit enter if reaction reference does not exist): ', 's');
disp(' ')

modName = input('Name of model reaction should be imported to: ', 's');
disp(' ')

fprintf('%s = addReaction(%s, ''%s'', %s, %s, %s, %s, %s, %s, ''%s'', '''', {}, %s, true);\n', modName, modName, rxnName, metString, Sstring, revFlag, lowerBound, upperBound, objCoeff, subSystem, metNameS);
fprintf('%s = addRxInfo(%s, %s, %s, %s, ''%s'', ''%s'', ''%s'', {});\n\n', modName, modName, metString, KEGGstring, metNameS, reactionName, rxRef, ECnum);
