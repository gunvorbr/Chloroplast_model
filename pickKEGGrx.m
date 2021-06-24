function [] = pickKEGGrx(model, rx, dir)

% Accepting KEGG reaction ID, and turning reaction into addReaction command
% that can be pasted on prompt to import reactions to model

%   by Gunvor Røkke, NTNU, 2021

disp(' ')
url = 'http://rest.kegg.jp/get/';

res = urlread(strcat(url, rx));
clear rx url

rxIDs = strfind(res, 'ENTRY') + 12; % Getting reaction ID
rxIDe = rxIDs + 5;
rxID = res(rxIDs:rxIDe);
rxRef = rxID;
clear rxIDs rxIDe
mode = input('Is reaction meant to be run in special organism mode? (1 = yes, 0 = no): ');
disp(' ')
if mode == 1
    org = input('Organism(s): ', 's');
    disp(' ')
    for i = 1:length(org)
        if i == 1
            orgStr = sprintf('@%s', org(i));
        else
            orgStr = sprintf('%s@%s', orgStr, org(i));
        end
    end
    rxID = sprintf('%s_%s', orgStr, rxID);
end
clear org i

KEGGIDss = strfind(res, 'EQUATION') + 12; % Collecting metabolite KEGG IDs
KEGGIDse = KEGGIDss;
while isempty(strfind(res(KEGGIDse), sprintf('\n')))
    KEGGIDse = KEGGIDse + 1;
end
KEGGIDstr = res(KEGGIDss:KEGGIDse);
clear KEGGIDss KEGGIDse

KEGGIDs = strsplit(KEGGIDstr);
clear KEGGIDstr

n = 1;

for i = 1:length(KEGGIDs) % Finding "trash" entries in KEGGIDs
    if ~isempty(strfind(char(KEGGIDs(i)), '+')) || isempty(char(KEGGIDs(i)))
        slettlist(n,1) = i;
        n = n + 1;
    end
end
for i = length(slettlist):-1:1 % Erasing "trash entries from KEGGIDs
    KEGGIDs(slettlist(i)) = [];
end
clear slettlist n

metSearch = findPattern(KEGGIDs, 'C'); % Working on S-matrix entry
Smat = zeros(1,size(metSearch,1));
n = 1;
i = 1;
arrow = findPattern(KEGGIDs, '>');
arrow = cell2mat(arrow(1,1));
while i <= length(KEGGIDs) % Filling S-matrix
    c = char(KEGGIDs(i));
    if ~isempty(strfind(char(KEGGIDs(i)), 'C'))
        if i < arrow
            Smat(1,n) = -1;
            n = n + 1;
        elseif i > arrow
            Smat(1,n) = 1;
            n = n + 1;
        end
    elseif isstrprop(c(1), 'digit')
        if i < arrow
            Smat(1,n) = -str2num(c);
            n = n + 1;
            i = i + 1;
        elseif i > arrow
            Smat(1,n) = str2num(c);
            n = n + 1;
            i = i + 1;
        end
    end
    i = i + 1;
end
clear c n i arrow KEGGIDs

KEGGmat = transpose(metSearch(:,2));

clear metSearch

if dir == -1 % Changing directionality of reaction
    revMat = -Smat;
    KEGGs = cell(1,length(KEGGmat));
    n = 1;
    for i = 1:length(revMat)
        if revMat(i) < 0
            Smat(n) = revMat(i);
            KEGGs(n) = KEGGmat(i);
            n = n + 1;
        end
    end
    for i = 1:length(revMat)
        if revMat(i) > 0
            Smat(n) = revMat(i);
            KEGGs(n) = KEGGmat(i);
            n = n + 1;
        end
    end
    KEGGmat = KEGGs;
    clear revMat KEGGs n
end

metIDs = cell(1,length(KEGGmat));
names = cell(1,length(KEGGmat));
comp = input('Choose compartment for reaction: ', 's');
disp(' ')

for i = 1:length(KEGGmat) % Filling metIDs and names
    if ~isempty(smatch(model.metKEGGID, char(KEGGmat(i)), 'exact')) % Metabolite is in model
        KEGGs = smatch(model.metKEGGID, char(KEGGmat(i)), 'exact');
        if length(KEGGs) > 1
            mets = model.mets(KEGGs);
            metRoots = cell(length(mets),1);
            for k = 1:length(mets)
                currentMet = char(mets(k));
                metRoots(k) = cellstr(currentMet(1:end-3));
            end
            clear mets currentMet k
            if length(unique(metRoots)) == 1
                oldMet = char(model.mets(KEGGs(1)));
            elseif length(unique(metRoots)) > 1
                fprintf('*** Metabolite %s ***\n\n', char(KEGGmat(i)))
                fprintf('Trying to assign metabolite ID to metabolite %s\n', char(KEGGmat(i)))
                fprintf('Multiple metabolites found in model under this KEGG ID. Please choose correct ID manually.\n')
                disp([transpose(num2cell(1:length(metRoots))), metRoots])
                indx = input('Index of correct metabolite (hit 0 if none of these are correct): ');
                disp(' ')
                if indx > 0
                    oldMet = char(model.mets(KEGGs(indx)));
                elseif indx == 0
                    ok = 0;
                    while ok == 0
                        pMetID = input('Propose a metabolite ID (compartment will be added automatically): ', 's');
                        pMetID = sprintf('%s[%s]', pMetID, comp);
                        if isempty(smatch(model.mets, pMetID, 'exact'))
                            oldMet = pMetID;
                            ok = 1;
                            disp('ID ok')
                            disp(' ')
                        else
                            fprintf('Metabolite ID is already found in model. Please try again\n')
                        end
                    end
                end
            end
        elseif length(KEGGs) == 1
            oldMet = char(model.mets(KEGGs));
        end
        metID = sprintf('%s[%s]', oldMet(1:end-3), comp);
        if exist('pMetID', 'var')
            oldName = input('Propose a metabolite name (compartment will be added automatically): ', 's');
            disp(' ')
            oldName = sprintf('%s [%s]', oldName, comp);
        else
            oldName = char(model.metNames(smatch(model.mets, oldMet, 'exact')));
        end
        clear pMetID
        if strcmp(oldName(end), ']') == 1
            name = sprintf('%s [%s]', oldName(1:end-4), comp);
        elseif strcmp(oldName(end-1), '_') == 0
            name = sprintf('%s [%s]', oldName(1:end-2, comp));
        end
        clear oldMet oldName
    elseif isempty(smatch(model.metKEGGID, char(KEGGmat(i)), 'exact')) % Metabolite is not in model
        fprintf('*** Metabolite %s ***\n\n', char(KEGGmat(i)))
        fprintf('Searching KEGG for metabolite. Name(s) found:\n')
        KEGGID = char(KEGGmat(i));
        metRes = urlread(strcat('http://rest.kegg.jp/get/', KEGGID));
        nameStart = strfind(metRes, 'NAME') + 12;
        nameEnd = nameStart;
        while ~strcmp(metRes(nameEnd), sprintf('\n')) || strcmp(metRes(nameEnd+1), ' ')
            nameEnd = nameEnd + 1;
        end
        nameEnd = nameEnd - 1;
        disp(' ')
        disp(metRes(nameStart:nameEnd))
        disp(' ')
        preName = input('Choose name for metabolite (compartment will be added automatically): ', 's');
        name = sprintf('%s [%s]', preName, comp);
        disp(' ')
        IDlock = 0;
        while IDlock == 0
            metID = input('Choose met ID for metabolite (compartment will be added automatically): ', 's');
            metID = sprintf('%s[%s]', metID, comp);
            if ~isempty(smatch(model.mets, metID, 'exact'))
                fprintf('Metabolite ID is already taken. Try again.\n\n')
            else
                IDlock = 1;
                disp('ID ok')
                disp(' ')
            end
        end
    end
    metIDs(i) = cellstr(metID);
    names(i) = cellstr(name);
end
clear dir i KEGGID KEGGs metID metRes name nameEnd nameStart oldMet oldName preName rootID




metList = metIDs; % Changing names on variables (because I'm lazy)
clear metIDs
metNames = names;
clear names
KEGGlist = KEGGmat;
clear KEGGmat
Slist = Smat;
clear Smat

if length(KEGGlist) > 1
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
elseif length(KEGGlist) == 1
    metString = sprintf('{''%s''}', char(metList(1)));
    metNameS = sprintf('{''%s''}', char(metNames(1)));
    KEGGstring = sprintf('{''%s''}', char(KEGGlist(1)));
    Sstring = sprintf('[%s]', num2str(Slist(1)));
end
clear j

disp(' ')
rev = input('To be, or not to be reversible? '); % Determining reversibility
disp(' ')

if (rev == 0) && (mode == 0) % Setting lower and upper bound
    lb = 0;
    ub = 1000;
elseif (rev == 1) && (mode == 0)
    lb = -1000;
    ub = 1000;
elseif (mode == 1) && ~strcmp(orgStr, '@N')
    lb = 0;
    ub = 0;
elseif (mode == 1) && strcmp(orgStr, '@N')
    ub = 1000;
    if rev == 0
        lb = 0;
    elseif rev == 1
        lb = -1000;
    end
end

dispSub = input('Sub system is next up. Display sub systems for inspiration? (yes = 1, no = 0): ');
disp(' ')
if dispSub == 1
    disp(char(unique(model.subSystems)))
    disp(' ')
end
subSystem = input('Sub-system? ', 's'); % Determining sub-system
disp(' ')

% Sjekk om NAME finnes før du starter løkke.
% Hvis ikke, display comment (hvis det finnes), og spør om å definere navn
% manuelt
if ~isempty(strfind(res, 'NAME'))
    fprintf('KEGG reaction name:\n')
    rxNameStart = strfind(res, 'NAME') + 12; % Setting reaction name
    rxNameSlutt = rxNameStart;
    while isempty(strfind(res(rxNameSlutt), sprintf('\n')))
        rxNameSlutt = rxNameSlutt + 1;
    end
    disp(' ')
    disp(res(rxNameStart:rxNameSlutt))
elseif isempty(strfind(res, 'NAME'))
    if ~isempty(strfind(res, 'COMMENT'))
        fprintf('Reaction name not present in KEGG. Reaction comment found:\n')
        comStart = strfind(res, 'COMMENT') + 12;
        comSlutt = comStart;
        while isempty(strfind(res(comSlutt), sprintf('\n')))
            comSlutt = comSlutt + 1;
        end
        disp(' ')
        disp(res(comStart:comSlutt))
    elseif isempty(strfind(res, 'COMMENT'))
        disp(' ')
        disp('No name information available')
    end
end
rxName = input('Set reaction name: ', 's');
disp(' ')
clear rxNameStart rxNameSlutt

% End edit

if strfind(res, 'ENZYME') % Setting EC number
    ECstart = strfind(res, 'ENZYME') + 12;
    ECslutt = length(res);
    ECstring = strsplit(res(ECstart:ECslutt));
    ECnum = char(ECstring(1));
    EC = 1;
else
    EC = 0;
end
clear ECstart ECslutt ECstring j

orgList = urlread('http://rest.kegg.jp/list/gn'); % yields list of organisms
orgName = input('Searching KEGG for genes. First name of organism you want to find genes for: ', 's');
disp(' ')
orgName = sprintf('%s%s', upper(orgName(1)), orgName(2:end)); % Making first letter uppercase (since strfind is case-sensitive)
orgHit = strfind(orgList, orgName);
if ~isempty(orgHit) % Organism exist in KEGG database
    startFlag = orgHit;
    while ~strcmp(orgList(startFlag:startFlag+3), 'gn:T')
        startFlag = startFlag - 1;
    end
    startFlag = startFlag + 10;
    endFlag = startFlag;
    while ~strcmp(orgList(endFlag), ',')
        endFlag = endFlag + 1;
    end
    endFlag = endFlag - 1;
    orgCode = upper(orgList(startFlag:endFlag));
    clear orgName orgList startFlag endFlag
    if EC == 1 % EC number exists for reaction
        EChit = urlread(strcat('http://rest.kegg.jp/get/', ECnum));
        geneHit = strfind(EChit, orgCode);
        if ~isempty(geneHit)
            geneStart = geneHit;
            while ~strcmp(EChit(geneStart), ':')
                geneStart = geneStart + 1;
            end
            geneStart = geneStart + 2;
            geneEnd = geneStart;
            while ~strcmp(EChit(geneEnd), sprintf('\n'))
                geneEnd = geneEnd + 1;
            end
            geneEnd = geneEnd - 1;
            geneStr = EChit(geneStart:geneEnd);
            genes = strsplit(geneStr);
            clear geneStr geneStart geneEnd orgCode EChit
            if length(genes) == 1 % Only one gene required for reaction to run
                fprintf('Gene for reaction found in KEGG: %s\n', char(genes))
                modify = input('Modify gene name? (yes = 1, no = 0): ');
                disp(' ')
                if modify == 1
                    genes = cellstr(input('New gene name: ', 's'));
                    disp(' ')
                end
                geneList = sprintf('{''%s''}', char(genes));
            elseif length(genes) > 1 % More than one gene required for reaction to run
                fprintf('Genes for reaction found in KEGG:\n')
                disp(transpose(genes))
                modify = input('Modify gene names? (yes = 1, no = 0): ');
                disp(' ')
                geneList = '{';
                for l = 1:length(genes)
                    if modify == 1
                        fprintf('Gene: %s\n', char(genes(l)))
                        genes(l) = cellstr(input('New gene name: ', 's'));
                        disp(' ')
                    end
                    geneList = sprintf('%s''%s'', ', geneList, char(genes(l)));
                end
                andOr = input('Are all genes necessary to run reaction (a), or only one (o)? ', 's');
                geneList = sprintf('%s''%s''}', geneList, andOr);
                clear andOr orgHit
            end
        end
    end
end

disp(' ')
disp(' ')

fprintf('chl = addReaction(chl, ''%s'', %s, %s, %i, %i, %i, 0, ''%s'', '''', {}, %s, true);\n', sprintf('%s_%s', rxID, comp), metString, Sstring, rev, lb, ub, subSystem, metNameS)
if EC == 1 && ~isempty(geneHit)
    fprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', ''%s'', ''%s'', %s);\n', metString, KEGGstring, metNameS, rxName, rxRef, ECnum, geneList)
elseif EC == 1 && isempty(geneHit)
    fprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', ''%s'', ''%s'', {});\n', metString, KEGGstring, metNameS, rxName, rxRef, ECnum)
elseif EC == 0 % Gene info also not existing
    fprintf('chl = addRxInfo(chl, %s, %s, %s, ''%s'', ''%s'', '''', {});\n', metString, KEGGstring, metNameS, rxName, rxRef)
end
disp(' ')
