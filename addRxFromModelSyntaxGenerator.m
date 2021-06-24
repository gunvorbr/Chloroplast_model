function syntaxRoll = addRxFromModelSyntaxGenerator(rxns, model, chl)

% Generates syntax for adding reactions to chloroplast model
% by picking reactions from other model (for example eyespot
% or organism-specific pathway).

% Input:
%   - List of reactions (#) to generate syntax for.
%   - Tag for indicating organism-specific reaction (optional).

% Output:
%   - Cell with syntax for adding reactions to chl model.

%   by Gunvor Røkke, NTNU, 2021

rxList = ReactionNames(model);
disp(' ')
org = input('Chlamy (C), Phaeo (P) or combination of organisms (CP)? ', 's');
disp(' ')
for i = 1:length(org)
    if i == 1
        orgStr = sprintf('@%s', org(i));
    else
        orgStr = sprintf('%s@%s', orgStr, org(i));
    end
end
org = orgStr;
clear orgStr
comp = input('Compartment symbol? ', 's');
syntaxRoll = cell(1,1);
n = 1;

fprintf('Overview of reactions to (potentially) be added:\n\n')
disp([model.rxns(rxns), model.rxnNames(rxns), rxList(rxns)])
fprintf('\n\n')

start = input('Are you ready to start? (1 = yes!): ');
disp(' ')
if start == 1
    load handel.mat
    sound(y,Fs)
end
clear y Fs
    

for a = 1:length(rxns)
    i = rxns(a);
    fprintf('\n\n*** %s ***\n\n', char(rxList(i)))
    addRx = input('Add reaction to chloroplast model? (1 = yes, 0 = no): ');
    disp(' ')
    if addRx == 1
        %% rxID og RxName
        fprintf('rxName: %s\n\n', char(model.rxnNames(i)))
        keep = input('Keep rxName (1) or add new (0)? ');
        disp(' ')
        if keep == 0
            rxName = input('New reaction name: ', 's');
            disp(' ')
        elseif keep == 1
            rxName = char(model.rxnNames(i));
        end
        rxID = char(model.rxns(i));
        fprintf('rxID: %s\n\n', char(model.rxns(i)))
        compInclude = input('Include compartment symbol in rxID? (yes = 1, no = 0): ');
        disp(' ')
        if compInclude == 1
            rxID = sprintf('%s_%s_%s', org, rxID, comp);
        elseif compInclude == 0
            rxID = sprintf('%s_%s', org, rxID);
        end
        % Checking that rxName and rxID is not already in chl
        if smatch(chl.rxnNames, rxName, 'exact')
            warning('Reaction name is already in model. Pick new reactionName')
            disp(' ')
            loop = 1;
            while loop == 1
                rxName = input('New reaction name: ', 's');
                if smatch(chl.rxns, rxName, 'exact')
                    warning('Reaction name still in model... Try again.')
                else
                    loop = 0;
                end
            end
        end
        if smatch(chl.rxns, rxID, 'exact')
            warning('Reaction ID is already in model. Pick new reaction ID')
            disp(' ')
            loop = 1;
            while loop == 1
                rxID = input('New reaction ID (organism / compartment will be added automatically): ', 's');
                if compInclude == 1
                    rxID = sprintf('%s_%s_%s', org, rxID, comp);
                elseif compInclude == 0
                    rxID = sprintf('%s_%s', org, rxID);
                end
                if smatch(chl.rxns, rxID, 'exact')
                    warning('Reaction name still in model... Try again.')
                else
                    loop = 0;
                end
            end
        end
        clear compInclude keep loop
        
        %% rxReference
        if isempty(model.rxnReferences{i})
            ref = 0;
        else
            ref = 1;
            rxRef = char(model.rxnReferences(i));
        end
        
        %% Sstring, metString, metNameS and KEGGstring
        
        indexMat = find(model.S(:,i) ~= 0);
        Smat = nonzeros(model.S(indexMat,i));
        metIDs = model.mets(indexMat);
        newMetIDs = cell(length(metIDs),1);
        metNames = model.metNames(indexMat);
        newMetNames = cell(length(metNames),1);
        KEGGIDs = model.metKEGGID(indexMat);
        
        % Ensuring that KEGG IDs are present (searching KEGG if it isn't
        missingKEGGs = cellfun(@isempty, KEGGIDs);
        for j = 1:length(KEGGIDs)
            if missingKEGGs(j) == 1
                metName = char(metNames(j));
                if strfind(metName(end), ']')
                    metName = metName(1:end-4);
                elseif strfind(metName(end-1), '_')
                    metName = metName(1:end-2);
                end
                fprintf('--- %s (%s) ---\n\n', metName, char(metIDs(j)));
                res = KmetChoose(metName);
                index = input('Index of correct KEGG-ID (search again / KEGG-ID not present = 0): ');
                disp(' ')
                if index == 0
                    again = input('Search again? (0 = no, 1 = yes): ');
                    disp(' ')
                    while again == 1
                        metName = input('New search term: ', 's');
                        res = KmetChoose(metName);
                        index = input('Index of correct KEGG-ID (search again / KEGG-ID not present = 0): ');
                        if index > 0
                            again = 0;
                        else
                            again = input('Search again = 1 / KEGG-ID not present = 0: ');
                        end
                        disp(' ')
                    end
                end
                if index > 0
                    index = num2str(index);
                    placeInTable = strmatch(index, res(:,1), 'exact');
                    if length(placeInTable) == 1
                        KEGGIDs(j) = res(placeInTable,2);
                    end
                elseif index == 0
                    search = input('Search model for metabolite? (yes = 1, no = 0): ');
                    if search == 1
                        fprintf('Metabolite name in model:    %s\n', metName)
                        fprintf('Metabolite ID in model:      %s\n\n', char(metIDs(j)))
                        while search == 1
                            array = input('Search in mets (1) or metNames (2)? ');
                            sTerm = input('Search term: ', 's');
                            if array == 1 % Searching model.mets
                                res = findPattern(chl.mets, sTerm);
                            elseif array == 2 % Searching model.metNames
                                res = findPattern(chl.metNames, sTerm);
                            end
                            found = input('Metabolite number (metabolite not found = 0): ');
                            if found > 0
                                KEGGIDs(j) = chl.metKEGGID(found);
                                search = 0;
                            elseif found == 0
                                search = input('Search again? (yes = 1, no = 0): ');
                            end
                        end
                    end
                end
            end
        end
        clear missingKEGGs j metName res index again placeInTable
        
        % Before translating to chl namespace, checking if reaction is exchange reaction
        origComps = cell(length(metIDs),1);
        for j = 1:length(metIDs)
            currID = char(metIDs(j));
            if strcmp(currID(end), ']')
                origComps(j) = cellstr(currID(end-1));
            elseif strcmp(currID(end-1), '_')
                origComps(j) = cellstr(currID(end));
            end
            clear currID
        end
        exchangeFlag = ~isequal(strcmp(origComps, char(origComps(1))), ones(length(origComps),1));
        if exchangeFlag == 1
            fprintf('\nReaction seems to be exchange reaction. Choose compartment symbols manually\n\n')
            chosenCompartments = cell(length(unique(origComps)),1);
            uniquomps = unique(origComps);
            for k = 1:length(uniquomps)
                fprintf('Choose compartment symbol for old compartment %s\n', char(uniquomps(k)))
                chosenCompartments(k) = cellstr(input('New compartment symbol: ', 's'));
            end
            disp(' ')
        end
        
        % Translating to format / language used in chl model
        for j = 1:length(metNames)
            currentKEGG = char(KEGGIDs(j));
            if ~isempty(currentKEGG)
                hits = smatch(chl.metKEGGID, currentKEGG, 'exact');
                if isempty(hits) % Metabolite is not already in chl model
                    if exchangeFlag == 0 % Creating ID / name for metabolite depending on exchangeFlag
                        currMetID = char(metIDs(j));
                        currMetName = char(metNames(j));
                        if strfind(currMetID(end), ']')
                            currMetID = sprintf('%s[%s]', currMetID(1:end-3), comp);
                        elseif strfind(currMetID(end-1), '_')
                            currMetID = sprintf('%s[%s]', currMetID(1:end-2), comp);
                        end
                        if strfind(currMetName(end-3:end-2), ' [')
                            currMetName = sprintf('%s [%s]', currMetName(1:end-4), comp);
                        else
                            warning('Help! (line 207)')
                        end
                    elseif exchangeFlag == 1
                        currMetID = char(metIDs(j));
                        currMetName = char(metNames(j));
                        compSymb = char(chosenCompartments(smatch(uniquomps, currMetID(end-1), 'exact')));
                        if strfind(currMetID(end), ']')
                            currMetID = sprintf('%s[%s]', currMetID(1:end-3), compSymb);
                        elseif strfind(currMetID(end-1), '_')
                            currMetID = sprintf('%s[%s]', currMetID(1:end-2), compSymb);
                        end
                        if strfind(currMetName(end-3:end-2), ' [')
                            currMetName = sprintf('%s [%s]', currMetName(1:end-4), comp);
                        else
                            warning('Help! (line 221)')
                        end
                    end
                    if smatch(chl.mets, currMetID, 'exact')
                        warning('Choose new ID for metabolite %s / %s', currMetID, currMetName)
                        compBackup = currMetID(end-1);
                        loop = 1;
                        while loop == 1
                            currMetID = input('Try new metabolite ID (compartment will be added automatically): ', 's');
                            currMetID = sprintf('%s[%s]', currMetID, compBackup);
                            if ~smatch(chl.mets, currMetID, 'exact')
                                loop = 0;
                            end
                        end
                    end
                    if smatch(chl.metNames, currMetName, 'exact')
                        warning('Choose new name for metabolite %s / %s', currMetID, currMetName)
                        compBackup = currMetName(end-1);
                        loop = 1;
                        while loop == 1
                            currMetName = input('Try new metabolite name (compartment will be added automatically): ', 's');
                            currMetName = sprintf('%s [%s]', currMetName, compBackup);
                            if ~smatch(chl.metNames, currMetName, 'exact')
                                loop = 0;
                            end
                        end
                    end
                    newMetNames(j) = cellstr(currMetName);
                    newMetIDs(j) = cellstr(currMetID);
                else % Metabolite is already in model
                    IDhits = chl.mets(hits);
                    if cellfun(@isempty, strfind(IDhits, ']')) ~= zeros(length(IDhits),1)
                        warning('Compartment symbols missing from metabolite IDs in reaction %s', rxName)
                    else
                        IDs = cell(length(IDhits),1);
                        for k = 1:length(IDs)
                            KID = char(IDhits(k));
                            IDs(k) = cellstr(KID(1:end-3));
                        end
                        if length(unique(IDs)) == 1
                            if exchangeFlag == 0
                                newMetID = sprintf('%s[%s]', char(unique(IDs)), comp);
                                newMetIDs(j) = cellstr(newMetID);
                            elseif exchangeFlag == 1
                                currMetID = char(metIDs(j));
                                compSymb = char(chosenCompartments(smatch(uniquomps, currMetID(end-1), 'exact')));
                                newMetID = sprintf('%s[%s]', char(unique(IDs)), compSymb);
                                newMetIDs(j) = cellstr(newMetID);
                            end 
                        end
                    end
                    nameHits = chl.metNames(hits);
                    if cellfun(@isempty, strfind(nameHits, ' [')) ~= zeros(length(nameHits),1)
                        warning('Compartment symbols missing from metabolite names in reaction %s', rxName)
                    else
                        names = cell(length(nameHits),1);
                        for k = 1:length(nameHits)
                            name = char(nameHits(k));
                            names(k) = cellstr(name(1:end-4));
                        end
                        if length(unique(names)) == 1
                            if exchangeFlag == 0
                                newMetName = sprintf('%s [%s]', char(unique(names)), comp);
                                newMetNames(j) = cellstr(newMetName);
                            elseif exchangeFlag == 1
                                currMetName = char(metNames(j));
                                compSymb = char(chosenCompartments(smatch(uniquomps, currMetName(end-1), 'exact')));
                                newMetName = sprintf('%s [%s]', char(unique(names)), compSymb);
                                newMetNames(j) = cellstr(newMetName);
                            end
                        end
                    end
                end
            elseif isempty(currentKEGG) % No KEGG ID present
                % Checking that old metabolite ID and metabolite name have compartments
                metID = char(metIDs(j));
                metName = char(metNames(j));
                if strfind(metID(end-1), '_') 
                    metID = sprintf('%s[%s]', metID(1:end-2), char(chosenCompartments(smatch(uniquomps, currMetID(end-1), 'exact'))));
                elseif strfind(metID(end-2), '[')
                    metID = sprintf('%s[%s]', metID(1:end-3), char(chosenCompartments(smatch(uniquomps, currMetID(end-1), 'exact'))));
                end
                if strfind(metName(end-3:end-2), ' [')
                    metName = sprintf('%s [%s]', metName(1:end-4), char(chosenCompartments(smatch(uniquomps, currMetID(end-1), 'exact'))));
                end
                % Checking that metID and metName are not already taken in chl. If one of them are, check if we are talking about the same metabolite...
                if smatch(chl.metNames, metName, 'exact')
                    nr = smatch(chl.metNames, metName, 'exact');
                    disp(' ')
                    dispMet(chl, nr)
                    disp(' ')
                    same = input('Is this the same metabolite as you are trying to add? (yes = 1, no = 0): ');
                    if same == 1
                        newMetIDs(j) = chl.mets(nr);
                        newMetNames(j) = chl.metNames(nr);
                        KEGGIDs(j) = chl.metKEGGID(nr);
                    end
                elseif smatch(chl.mets, metID, 'exact')
                    nr = smatch(chl.mets, metID, 'exact');
                    disp(' ')
                    dispMet(chl, nr)
                    disp(' ')
                    same = input('Is this the same metabolit as you are trying to add? (yes = 1, no = 0): ');
                    if same == 1
                        newMetIDs(j) = chl.mets(nr);
                        newMetNames(j) = chl.metNames(nr);
                        KEGGIDs(j) = chl.metKEGGID(nr);
                    end
                end
                if isempty(newMetNames{j})
                    if isempty(smatch(chl.mets, metID, 'exact')) || isempty(smatch(chl.metNames, metName, 'exact'))
                        newMetIDs(j) = metIDs(j);
                        newMetNames(j) = metNames(j);
                        KEGGIDs(j) = cellstr('');
                    else
                        newMetIDs(j) = cellstr(input('Choose metabolite ID: ', 's'));
                        newMetNames(j) = cellstr(input('Choose metabolite name: ', 's'));
                        KEGGIDs(j) = cellstr('');
                        search = 1;
                        while search == 1
                            if isempty(smatch(chl.mets, char(newMetID(j)), 'exact')) || isempty(smatch(chl.metNames, char(newMetNames(j)), 'exact'))
                                search = 0;
                            else
                                warning('Name or ID taken! Please choose again!\n\n')
                                newMetIDs(j) = cellstr(input('Choose metabolite ID: ', 's'));
                                newMetNames(j) = cellstr(input('Choose metabolite name: ', 's'));
                            end
                        end
                    end
                end
            end
        end
        disp(' ')
        % Checking
        disp('Old vs. new met IDs:')
        disp([metIDs, newMetIDs])
        disp(' ')
        disp('Old vs. new met names:')
        disp([metNames, newMetNames])
        disp(' ')
        %% Generating Sstring, metString and metNames
        if length(Smat) == 1
            Sstring = sprintf('%i', Smat);
            metString = sprintf('{''%s''}', char(newMetIDs));
            metNameS = sprintf('{''%s''}', char(newMetNames));
            KEGGstring = sprintf('{''%s''}', char(KEGGIDs));
        else
            Sstring = sprintf('[%i', Smat(1));
            metString = sprintf('{''%s''', char(newMetIDs(1)));
            metNameS = sprintf('{''%s''', char(newMetNames(1)));
            KEGGstring = sprintf('{''%s''', char(KEGGIDs(1)));
            for j = 2:length(Smat)
                Sstring = sprintf('%s, %i', Sstring, Smat(j));
                metString = sprintf('%s, ''%s''', metString, char(newMetIDs(j)));
                metNameS = sprintf('%s, ''%s''', metNameS, char(newMetNames(j)));
                KEGGstring = sprintf('%s, ''%s''', KEGGstring, char(KEGGIDs(j)));
            end
            Sstring = sprintf('%s]', Sstring);
            metString = sprintf('%s}', metString);
            metNameS = sprintf('%s}', metNameS);
            KEGGstring = sprintf('%s}', KEGGstring);
        end
        %% Reversity, lower bound, upper bound
        rev = model.rev(i);
        if ~strcmp(org, '@N')
            lb = 0;
            ub = 0;
        elseif strcmp(org, '@N')
            lb = model.lb(i);
            ub = model.ub(i);
        end
        if strcmp(org, '@N') && (((rev == 1) && (lb == 0)) || ((rev == 0) && (lb < 0)))
            warning('Muffens is going on with reversibility and lower bound :-(')
        end
        %% Sub system
        fprintf('Current sub system: %s\n', char(model.subSystems(i)))
        disp(' ')
        display = input('Display chl sub system for inspiration? (0 = no, 1 = yes): ');
        if display == 1
            disp(unique(chl.subSystems))
            disp(' ')
        end
        keep = input('Keep sub system (1) or add new (0): ');
        if keep == 1
            subSystem = char(model.subSystems(i));
        elseif keep == 0
            subSystem = input('Choose sub system: ', 's');
        end
        %% EC number
        if isempty(model.rxnECNumbers{i})
            EC = 0;
        else
            ECnum = char(model.rxnECNumbers(i));
            EC = 1;
        end
        %% Genes
        clear geneIndx genes geneNames andOr geneList
        geneIndx = find(model.rxnGeneMat(i,:) ~= 0);
        if ~isempty(geneIndx)
            genes = 1;
            geneNames = model.genes(geneIndx);
            if ~isempty(strfind(lower(char(model.grRules(i))), 'and')) && isempty(strfind(lower(char(model.grRules(i))), 'or'))
                andOr = 'a';
            elseif ~isempty(strfind(lower(char(model.grRules(i))), 'or')) && isempty(strfind(lower(char(model.grRules(i))), 'and'))
                andOr = 'o';
            elseif ~isempty(strfind(lower(char(model.grRules(i))), 'and')) && ~isempty(strfind(lower(char(model.grRules(i))), 'or'))
                andOr = 'a';
            end
            geneList = '{';
            for l = 1:length(geneNames)
                geneList = sprintf('%s''%s'', ', geneList, char(geneNames(l)));
            end
            geneList = sprintf('%s''%s''}', geneList, andOr);
        elseif isempty(geneIndx)
            genes = 0;
        end
        %% Generating syntax and adding to syntax roll
        syntaxRoll(n,1) = cellstr(sprintf('chl = addReaction(chl, ''%s'', %s, %s, %i, %i, %i, 0, ''%s'', '''', {}, %s, true);', rxID, metString, Sstring, rev, lb, ub, subSystem, metNameS));
        n = n + 1;
        if (EC == 0) && (ref == 0) && (genes == 0)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', '''', '''', {});', metString, KEGGstring, metNameS, rxName));
            n = n + 1;
        elseif (EC == 0) && (ref == 0) && (genes == 1)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', '''', '''', %s);', metString, KEGGstring, metNameS, rxName, geneList));
            n = n + 1;
        elseif (EC == 0) && (ref == 1) && (genes == 0)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', ''%s'', '''', {});', metString, KEGGstring, metNameS, rxName, rxRef));
            n = n + 1;
        elseif (EC == 1) && (ref == 0) && (genes == 0)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', '''', ''%s'', {});', metString, KEGGstring, metNameS, rxName, ECnum));
            n = n + 1;
        elseif (EC == 1) && (ref == 1) && (genes == 0)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', ''%s'', ''%s'', {});', metString, KEGGstring, metNameS, rxName, rxRef, ECnum));
            n = n + 1;
        elseif (EC == 1) && (ref == 0) && (genes == 1)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', '''', ''%s'', %s);', metString, KEGGstring, metNameS, rxName, ECnum, geneList));
            n = n + 1;
        elseif (EC == 0) && (ref == 1) && (genes == 1)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', ''%s'', '''', %s);', metString, KEGGstring, metNameS, rxName, rxRef, geneList));
            n = n + 1;
        elseif (EC == 1) && (ref == 1) && (genes == 1)
            syntaxRoll(n,1) = cellstr(sprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', ''%s'', ''%s'', %s);', metString, KEGGstring, metNameS, rxName, rxRef, ECnum, geneList));
            n = n + 1;
        end
        disp(' ')
        fprintf('Syntax roll so far (after rx %i):\n', i);
        disp(' ')
        disp(char(syntaxRoll))
        disp(' ')
        syntaxRoll(n,1) = cellstr(' ');
        n = n + 1;
    end
end
