function [model, missingMets] = plugAndPlay(cyt, chl)

% Merges model without chloroplast compartment with chloroplast model

% Additional scripts required to run plugAndPlay.m:
%   - metKEGGIDsearch.m
%   - fillKEGGIDholes.m
%   - KEGGIDfirstAid.m
%   - smatch.m
%   - ReactionNames.m

%   by Gunvor Røkke, NTNU, 2021

%% 0 - Saving backup of models

fprintf('\nThis script will make changes to input models\n\n')
backup = input('Do you want to save a backup version of the models? (yes = 1, no = 0): ');
disp(' ')
if backup == 1
    cyt_filename = sprintf('%s_Backup_cytosol_model.mat', date);
    chl_filename = sprintf('%s_Backup_chloroplast_model.mat', date);
    save(cyt_filename, 'cyt')
    save(chl_filename, 'chl')
end
clear cyt_filename chl_filename backup

%% 1 - Deleting biomass reactions from chl model

bioList = smatch(chl.rxnNames, 'biomass');
bioNames = chl.rxns(bioList);
for i = 1:length(bioNames)
    chl = removeRxns(chl, bioNames(i));
end
clear bioList bioNames

%% 2 - Checking compartment symbol overlap

% Checking that the compartment symbols of the two models to be merged do
% not overlap

met = char(cyt.mets(1)); % Determining where compartment symbol is and changing format if necessary
if strcmp(met(end-1), '_')
    for i = 1:length(cyt.mets) % Changing format of cyt.mets vector
        oMet = char(cyt.mets(i));
        comp = oMet(end);
        nMet = sprintf('%s[%s]', oMet(1:end-2), comp);
        cyt.mets(i) = cellstr(nMet);
    end
    clear oMet comp nMet
elseif ~strcmp(met(end), ']')
    fprintf('Example of metabolite ID: %s\n\n', met)
    compSymb = input('Compartment symbol in example above: ', 's');
    disp(' ')
    fullID = input('Full metabolite ID in example above: ', 's');
    disp(' ')
    pos = input('Compartment symbol situated at the beginning (1) or the end (2)? ');
    disp(' ')
    if pos == 1
        x = 1;
        while ~strcmp(met(x), compSymb)
            x = x + 1;
        end
        IDstart = strfind(met, fullID);
    elseif pos == 2
        x = length(met);
        while ~strcmp(met(x), compSymb)
            x = x - 1;
        end
        IDend = strfind(met, fullID) + length(fullID) - 1;
    end
    clear compSymb
    for i = 1:length(cyt.mets) % Changing format of cyt.mets vector
        oMet = char(cyt.mets(i));
        comp = oMet(x);
        if pos == 1
            nMet = sprintf('%s[%s]', oMet(IDstart:end), comp);
        elseif pos == 2
            nMet = sprintf('%s[%s]', oMet(1:IDend), comp);
        end
        cyt.mets(i) = cellstr(nMet);
    end
    clear oMet comp nMet IDstart IDend
end
clear met pos fullID

% Changing format of compartment symbols in cytosol model if compFlag ~= 1

metNameCheck = char(cyt.metNames(1));
if (strcmp(metNameCheck(end), ']') == 0) || (strcmp(metNameCheck(end-2), '_') > 0) || (strcmp(metNameCheck(end-1), '_') > 0)
    if strcmp(metNameCheck(end-2), '_')
        nameFlag = 3;
    elseif strcmp(metNameCheck(end-1), '_')
        nameFlag = 2;
    else
        nameFlag = 1;
    end
    for i = 1:length(cyt.metNames)
        name = char(cyt.metNames(i));
        met = char(cyt.mets(i));
        comp = met(end-1);
        if nameFlag == 1
            nName = sprintf('%s [%s]', name, comp);
        elseif nameFlag == 2
            nName = sprintf('%s [%s]', name(1:end-2), comp);
        elseif nameFlag == 3
            nName = sprintf('%s [%s]', name(1:end-3), comp);
        end
        cyt.metNames(i) = cellstr(nName);
    end
end
clear metNameCheck nameFlag name met comp nName

% Checking that [c] is the only overlapping compartment between the models

cytComps = cell.empty(0,1); % Making empty vector for compartments of cytosol model metabolites
for i = 1:length(cyt.mets) % Filling cytosol model metabolite vector
    met = char(cyt.mets(i));
    if strcmp(met(end), ']')
        cytComps(i,1) = cellstr(met(end-1));
    else
        warning('Something must have gone wrong with conversion of IDs to correct format...')
    end
end
clear met

cytComps = unique(cytComps); % Reducing sice of cytosol model metabolite vector

if ~isequal(strcmp(cytComps, 'h'), zeros(length(cytComps),1)) % Checking if chloroplast symbol h is found amongst compartment symbols of cytosol model
    rxList = ReactionNames(cyt);
    hRxns = cell.empty(0,1);
    n = 1;
    for j = 1:length(rxList)
        currentEq = char(rxList(j));
        if contains(currentEq, '[h]')
            hRxns(n,1) = cyt.rxns(j);
            n = n + 1;
        end
    end
    bioRx = 1;
    for j = 1:length(hRxns)
        if ~contains(lower(char(hRxns(j))), 'bio')
            bioRx = 0;
        end
    end
    if bioRx == 0
        fprintf('Symbol for chloroplast in chloroplast model (h) is also found in cytosol model\n')
        newChloroplastSymb = input('Choose new chloroplast symbol: ', 's');
        disp(' ')
        for i = 1:length(chl.mets) % Changing chloroplast symbol in chloroplast model
            met = char(chl.mets(i));
            if strcmp(met(end-1), 'h')
                newMetID = sprintf('%s[%s]', met(1:end-3), newChloroplastSymb);
                chl.mets(i) = cellstr(newMetID);
                oldMetName = char(chl.metNames(i));
                newMetName = sprintf('%s [%s]', oldMetName(1:end-4), newChloroplastSymb);
                chl.metNames(i) = cellstr(newMetName);
            end
        end
    end
    clear newChloroplastSymb met newMetID oldMetName newMetName currentEq hRxns bioRx n j
end
if ~isequal(strcmp(cytComps, 'u'), zeros(length(cytComps),1)) % Checking if thylakoid symbol u is found amongst compartment symbols of cytosol model
    rxList = ReactionNames(cyt);
    uRxns = cell.empty(0,1);
    n = 1;
    for j = 1:length(rxList)
        currentEq = char(rxList(j));
        if contains(currentEq, '[u]')
            uRxns(n,1) = cyt.rxns(j);
            n = n + 1;
        end
    end
    bioRx = 1;
    for j = 1:length(uRxns)
        if ~contains(lower(char(uRxns(j))), 'bio')
            bioRx = 0;
        end
    end
    if bioRx == 0
        fprintf('Symbol for thylakoid in chloroplast model (u) is also found in cytosol model\n')
        newThylakoidSymb = input('Choose new thylakoid symbol: ', 's');
        disp(' ')
        for i = 1:length(chl.mets) % Changing thylakoid symbol n chloroplast model
            met = char(chl.mets(i));
            if strcmp(met(end-1), 'u')
                newMetID = sprintf('%s[%s]', met(1:end-3), newThylakoidSymb);
                chl.mets(i) = cellstr(newMetID);
                oldMetName = char(chl.metNames(i));
                newMetName = sprintf('%s [%s]', oldMetName(1:end-4), newThylakoidSymb);
                chl.metNames(i) = cellstr(newMetName);
            end
        end
    end
    clear newThylakoidSymb met newMetID oldMetName newMetName currentEq uRxns bioRx j n
end
if ~isequal(strcmp(cytComps, 's'), zeros(length(cytComps),1)) % Checking if eyespot symbol s is found amongst compartment symbols of cytosol model
    rxList = ReactionNames(cyt);
    sRxns = cell.empty(0,1);
    n = 1;
    for j = 1:length(rxList)
        currentEq = char(rxList(j));
        if contains(currentEq, '[s]')
            sRxns(n,1) = cyt.rxns(j);
            n = n + 1;
        end
    end
    bioRx = 1;
    for j = 1:length(sRxns)
        if ~contains(lower(char(sRxns(j))), 'bio')
            bioRx = 0;
        end
    end
    if bioRx == 0
        fprintf('Symbol for eyespot in chloroplast model (s) is also found in cytosol model\n')
        newEyespotSymb = input('Choose new eyespot symbol: ', 's');
        disp(' ')
        for i = 1:length(chl.mets) % Changing eyespot symbol n chloroplast model
            met = char(chl.mets(i));
            if strcmp(met(end-1), 's')
                newMetID = sprintf('%s[%s]', met(1:end-3), newEyespotSymb);
                chl.mets(i) = cellstr(newMetID);
                oldMetName = char(chl.metNames(i));
                newMetName = sprintf('%s [%s]', oldMetName(1:end-4), newEyespotSymb);
                chl.metNames(i) = cellstr(newMetName);
            end
        end
    end
    clear newThylakoidSymb met newMetID oldMetName newMetName currentEq sRxns bioRx n j
end
clear cytComps met i

%% 3 - Making sure cytosol model contains KEGG IDs for metabolites

KEGGsPresent = ~cellfun(@isempty, cyt.metKEGGID); % Checking how many metabolites in cytosol model that has a KEGG ID
KEGGID = length(find(KEGGsPresent == 1)); % Determining amount of metabolites that do or do not have KEGG IDs
noKEGGID = length(find(KEGGsPresent == 0));

if (KEGGID / (KEGGID + noKEGGID)) < 0.6
    fprintf('KEGG IDs of metabolites are necessary to connect similar metabolites of the two models\n');
    fprintf('%i out of %i metabolites in the cytosol model are not associated with a KEGG ID.\n\n', noKEGGID, length(cyt.mets));
    fill = input('Fill in KEGG-IDs? (yes (use detailed method) = 1, yes (use quick method) = 2, no = 0): ');
    disp(' ')
    if fill == 1 % Detailed fill-mode
        cyt = metKEGGIDsearch(cyt); % Attemting to find more metabolite KEGG IDs
        
        stillEmptyKEGGs = 0; % Checking if there are more KEGG ID holes to be filled
        for i = 1:length(cyt.metKEGGID)
            if isempty(cyt.metKEGGID{i})
                stillEmptyKEGGs = 1;
            end
        end
        if stillEmptyKEGGs == 1 % Filling more KEGG IDs
            fprintf('Attempting to fill in the last missing KEGG IDs. Hold on tight...\n\n')
            cyt = fillKEGGIDholes(cyt);
        end
        
        stillEmptyKEGGs = 0; % Checking if there are metabolites without KEGG ID after two rounds of filling
        for i = 1:length(cyt.metKEGGID)
            if isempty(cyt.metKEGGID{i})
                stillEmptyKEGGs = 1;
            end
        end
        
        if stillEmptyKEGGs == 0 % Celebrating :-D
            fprintf('Every metabolite now has a KEGG ID! Hooray :-D\n\n')
            load handel.mat
            sound(y,Fs)
            clear y Fs
        end
    elseif fill == 2
        fprintf('Attempting to find KEGG-IDs only for the metabolites being transported between the chloroplast and the cytoplasm\n\n')
        cyt = KEGGIDfirstAid(cyt, chl);
    end
    if fill > 0
        backup = input('Save backup of cytosol model now that KEGG IDs have been added? (yes = 1, no = 0): ');
        disp(' ')
        if backup == 1
            filename = input('Filename? (.mat will be added automatically): ', 's');
            disp(' ')
            filename = sprintf('%s.mat', filename);
            save(filename, 'cyt');
        end
    end
end

clear i fill KEGGID KEGGsPresent noKEGGID stillEmptyKEGGs backup

%% 4 - Translating exchange metabolite of chloroplast model to cytosol namespace

fprintf('Translating exchange metabolites in chloroplast model into cytosol model namespace\n\n')

cMets = double.empty(0,1); % Generating list of cytosol metabolites in chloroplast model
n = 1;
for i = 1:length(chl.mets)
    met = char(chl.mets(i));
    if strcmp(met(end-1), 'c')
        cMets(n,1) = i;
        n = n + 1;
    end
end
pMets = zeros(length(cMets),1); % Identifying plastid version of cytosol metabolites
for a = 1:length(cMets)
    i = cMets(a);
    met = char(chl.mets(i));
    pMet = sprintf('%s[h]', met(1:end-3));
    if ~isempty(smatch(chl.mets, pMet, 'exact'))
        if length(smatch(chl.mets, pMet, 'exact')) == 1
            pMets(a,1) = smatch(chl.mets, pMet, 'exact');
        end
    end
end
uMets = zeros(length(cMets),1); % Determining if exchange metabolite is also found in thylakoid
for a = 1:length(cMets)
    i = cMets(a);
    met = char(chl.mets(i));
    uMet = sprintf('%s[u]', met(1:end-3));
    if ~isempty(smatch(chl.mets, uMet, 'exact'))
        if length(smatch(chl.mets, uMet, 'exact')) == 1
            uMets(a,1) = smatch(chl.mets, uMet, 'exact');
        end
    end
end
sMets = zeros(length(cMets),1); % Determining if exchange metabolite is also found in eyespot
for a = 1:length(cMets)
    i = cMets(a);
    met = char(chl.mets(i));
    sMet = sprintf('%s[s]', met(1:end-3));
    if ~isempty(smatch(chl.mets, sMet, 'exact'))
        if length(smatch(chl.mets, sMet, 'exact')) == 1
            sMets(a,1) = smatch(chl.mets, sMet, 'exact');
        end
    end
end

for i = length(pMets):-1:1 % Removing cytosol metabolites that does not have a plastid counterpart
    if pMets(i) == 0;
        pMets(i) = [];
        cMets(i) = [];
        sMets(i) = [];
        uMets(i) = [];
    end
end
clear i met n pMet sMet uMet

missingMets = cell.empty(0,2); % Making cell array for chl model c-mets not present in cyt model
n = 1;

for a = 1:length(cMets) % Searching for cMets in cytosol model, translating cMets and pMets (and if available; sMets) into cytosol namespace
    i = cMets(a);
    k = pMets(a);
    u = uMets(a);
    s = sMets(a);
    metName = char(chl.metNames(i));
    KEGGID = char(chl.metKEGGID(i));
    if ~contains(lower(metName), 'chloroplast') % Regular metabolite, not biomass metabolite
        chlMet = char(chl.mets(i));
        metMatch = smatch(cyt.metKEGGID, KEGGID, 'exact');
        if ~isempty(metMatch)
            constName = 1;
            for j = 1:length(metMatch) % Going through metabolite matches based on KEGG ID from cyt model. If metabolite roots are different, the metabolite in question is probabaly a fatty acid. If that is the case, program lets user pick correct metabolite.
                if j == 1
                    refMet = char(cyt.mets(metMatch(j)));
                    refRoot = refMet(1:end-3);
                else
                    checkMet = char(cyt.mets(metMatch(j)));
                    checkRoot = checkMet(1:end-3);
                    if ~strcmp(refRoot, checkRoot)
                        constName = 0;
                    end
                end
            end
            if constName == 0 % Cyt metabolite IDs are on different format. Choosing manually
                fprintf('* Altering namespace for chloroplast metabolite %s\n\n', chlMet)
                fprintf('Metabolite hits in cytosol model:\n')
                disp([transpose(num2cell(1:length(metMatch))), cyt.mets(metMatch)])
                correctMet = input('Index of correct metabolite (if correct metabolite is not present, press 0): ');
                disp(' ')
                if correctMet > 0
                    cytCmet = metMatch(correctMet);
                elseif correctMet == 0
                    cytCmet = 0;
                end
            elseif constName == 1 % Cyt metabolite IDs are on the same format
                if length(metMatch) == 1 % Searching for cytosol metabolite in cytosol model
                    cytMet = char(cyt.mets(metMatch));
                    if strcmp(cytMet(end-1), 'c')
                        cytCmet = metMatch;
                    else
                        cytCmet = 0;
                    end
                else
                    cytCmet = double.empty(0,1);
                    m = 1;
                    for j = 1:length(metMatch)
                        cytMet = char(cyt.mets(metMatch(j)));
                        if strcmp(cytMet(end-1), 'c')
                            cytCmet(m,1) = metMatch(j);
                            m = m + 1;
                        end
                    end
                    if length(cytCmet) > 1
                        fprintf('Input required for chloroplast metabolite %s\n\n', chlMet)
                        fprintf('Choose correct cytosol metabolite:\n')
                        disp([transpose(num2cell(1:length(cytCmet))), cyt.mets(cytCmet)])
                        cytCmet = input('Choose index of correct metabolite: ');
                        disp(' ')
                    end
                end
            end
            if cytCmet > 0 % Changing name of cytosol metabolite and corresponding plastid metabolite in chloroplast model
                % Changing cytosol metabolite ID
                oldCid = char(chl.mets(i));
                if ~strcmp(cyt.mets(cytCmet), chl.mets(i)) && ~isempty(smatch(chl.mets, char(cyt.mets(cytCmet)), 'exact'))
                    fprintf('* Metabolite %i of %i: Cannot change name of metabolite %s to %s. %s found in chl model already\n', a, length(cMets), oldCid, char(cyt.mets(cytCmet)), char(cyt.mets(cytCmet)))
                    changeOK = 0;
                else
                    chl.mets(i) = cyt.mets(cytCmet);
                    fprintf('* Metabolite %i of %i: Metabolite ID changed from %s to %s\n', a, length(cMets), oldCid, char(chl.mets(i)))
                    newCid = char(chl.mets(i));
                    changeOK = 1;
                end
                % Changing plastid metabolite ID
                oldPid = char(chl.mets(k));
                newPid = sprintf('%s[h]', newCid(1:end-3));
                if changeOK == 1
                    chl.mets(k) = cellstr(newPid);
                    fprintf('  Metabolite %i of %i: Metabolite ID changed from %s to %s\n', a, length(cMets), oldPid, char(chl.mets(k)))
                else
                    fprintf('* Metabolite %i of %i: Cannot change name of metabolite %s to %s. %s found in chl model already\n', a, length(cMets), oldCid, char(cyt.mets(cytCmet)), char(cyt.mets(cytCmet)))
                end
                % Changing thylakoid metabolite ID
                if u ~= 0
                    oldTid = char(chl.mets(u));
                    newTid = sprintf('%s[u]', newCid(1:end-3));
                    if changeOK == 1
                        chl.mets(u) = cellstr(newTid);
                        fprintf('  Metabolite %i of %i: Metabolite ID changed from %s to %s\n', a, length(cMets), oldTid, char(chl.mets(u)))
                    else
                        fprintf('* Metabolite %i of %i: Cannot change name of metabolite %s to %s. %s found in chl model already\n', a, length(cMets), oldTid, newTid, newTid)
                    end
                end
                % Changing eyespot metabolite ID
                if s ~= 0
                    oldSid = char(chl.mets(s));
                    newSid = sprintf('%s[s]', newCid(1:end-3));
                    if changeOK == 1
                        chl.mets(s) = cellstr(newSid);
                        fprintf('  Metabolite %i of %i: Metabolite ID changed from %s to %s\n', a, length(cMets), oldSid, char(chl.mets(s)))
                    else
                        fprintf('* Metabolite %i of %i: Cannot change name of metabolite %s to %s. %s found in chl model already\n', a, length(cMets), oldSid, newSid, newSid)
                    end
                end
                % Changing cytosol metabolite name
                oldCname = char(chl.metNames(i));
                chl.metNames(i) = cyt.metNames(cytCmet);
                fprintf('  Metabolite %i of %i: Metabolite name changed from %s to %s\n', a, length(cMets), oldCname, char(chl.metNames(i)))
                newCname = char(chl.metNames(i));
                % Changing plastid metabolite name
                oldPname = char(chl.metNames(k));
                newPname = sprintf('%s [h]', newCname(1:end-4));
                chl.metNames(k) = cellstr(newPname);
                fprintf('  Metabolite %i of %i: Metabolite name changed from %s to %s\n', a, length(cMets), oldPname, char(chl.metNames(k)))
                if u ~= 0 % Changing thylakoid metabolite name
                    oldTname = char(chl.metNames(u));
                    newTname = sprintf('%s [u]', newCname(1:end-4));
                    chl.metNames(u) = cellstr(newTname);
                    fprintf('  Metabolite %i of %i: Metabolite name changed from %s to %s\n', a, length(cMets), oldTname, char(chl.metNames(u)))
                end
                if s ~= 0 % Changing eyespot metabolite name
                    oldSname = char(chl.metNames(s));
                    newSname = sprintf('%s [s]', newCname(1:end-4));
                    chl.metNames(s) = cellstr(newSname);
                    fprintf('  Metabolite %i of %i: Metabolite name changed from %s to %s\n', a, length(cMets), oldSname, char(chl.metNames(s)))
                end
                disp(' ')
                clear oldCid newCid oldPid newPid oldTid newTid oldSid newSid oldCname newCname oldPname newPname oldTname newTname oldSname newSname
            elseif cytCmet == 0 % Metabolite is not present in cytosol compartment of cytosol model
                missingMets(n,1) = cellstr(metName);
                missingMets(n,2) = cellstr(KEGGID);
                n = n + 1;
                fprintf('* Metabolite %i of %i:\n', a, length(cMets))
                fprintf('  Metabolite %s is present in cytosol model, but does not seem to be present in the cytosol compartment\n', metName(1:end-4))
                fprintf('  Metabolite has been added to list of missing metabolites\n\n')
            end
        elseif isempty(metMatch) % Metabolite might not be present in cytosol model
            missingMets(n,1) = cellstr(metName);
            missingMets(n,2) = cellstr(KEGGID);
            n = n + 1;
            fprintf('* Metabolite %i of %i:\n', a, length(cMets))
            fprintf('  Metabolite %s is seemingly not present in cytosol model\n', metName(1:end-4))
            fprintf('  Metabolite has been added to list of missing metabolites\n\n')
        end
    end
end
clear cytCmet cytMet a i j k m n s u KEGGID chlMet metMatch metName constName refMet refRoot checkMet checkRoot correctMet changeOK cMets pMets uMets sMets

%% 5 - Merging models

fprintf('Merging models. Please wait...\n\n')
model = cyt;

for j = 1:length(chl.rxns) % Looping through chloroplast model reactions
    n = length(model.rxns) + 1; % Index of next reaction from chl model to be added to model
    % rxns
    model.rxns(n) = chl.rxns(j);
    % rxnName
    model.rxnNames(n) = chl.rxnNames(j);
    % subSystems
    model.subSystems(n) = chl.subSystems(j);
    % rxnECNumbers
%    model.rxnECNumbers(n) = chl.rxnECNumbers(j);
    % rev
    model.rev(n) = chl.rev(j);
    % ub
    model.ub(n) = chl.ub(j);
    % lb
    model.lb(n) = chl.lb(j);
    % c
    model.c(n) = 0;
    % confidenceScores
%    model.confidenceScores(n) = chl.confidenceScores(j);
    % rules
    model.rules(n) = chl.rules(j);
    % grRules
    model.grRules(n) = chl.grRules(j);
    % rxnReferences
%    model.rxnReferences(n) = chl.rxnReferences(j);
    % rxnNotes
%    model.rxnNotes(n) = chl.rxnNotes(j);
    
    % Handling metabolites involved in reaction
    metMatrix = [find(chl.S(:,j) ~= 0), nonzeros(chl.S(:,j)), zeros(length(nonzeros(chl.S(:,j))),1), zeros(length(nonzeros(chl.S(:,j))),1)]; % met index in chl, S-coeff in chl, met index in model (for now set to zero), metabolite found in model already?
    metIDs = chl.mets(metMatrix(:,1));
    m = length(model.mets) + 1;
    for i = 1:length(metIDs)
        met = char(metIDs(i));
        modelMetIndex = smatch(model.mets, met, 'exact');
        if ~isempty(modelMetIndex)
            if length(modelMetIndex) == 1
                metMatrix(i,3) = modelMetIndex;
                metMatrix(i,4) = 1;
            else
                warning('Metabolite %s is found more than once in metabolite vector of model', met)
                disp(' ')
            end
        elseif isempty(modelMetIndex)
            metMatrix(i,3) = m;
            m = m + 1;
        end
    end
    
    for i = 1:size(metMatrix,1)
        h = metMatrix(i,1); % Metabolite index in chl model
        c = metMatrix(i,3); % Metabolite index in model
        if metMatrix(i,4) == 0
            % mets
            model.mets(c) = chl.mets(h);
            % metNames
            model.metNames(c) = chl.metNames(h);
            % metKEGGID
%            model.metKEGGID(c) = chl.metKEGGID(h);
            % b
%            model.b(c) = 0;
            % metFormulas
%            model.metFormulas(c) = chl.metFormulas(h);
            % metChEBIID
%            model.metChEBIID(c) = chl.metChEBIID(h);
            % metPubChemID
%            model.metPubChemID(c) = chl.metPubChemID(h);
            % metInChIString
%            model.metInChIString(c) = chl.metInChIString(h);
            % metCharge
%            model.metCharge(c) = chl.metCharge(h);
        elseif metMatrix(i,4) == 1
%            if isempty(model.metKEGGID{c})
                % metKEGGID
%                model.metKEGGID(c) = chl.metKEGGID(h);
%            end
%            if isempty(model.metFormulas{c})
                % metFormulas
%                model.metFormulas(c) = chl.metFormulas(h);
%            end
%            if isempty(model.metChEBIID{c})
                % metChEBIID
%                model.metChEBIID(c) = chl.metChEBIID(h);
%            end
%            if isempty(model.metPubChemID{c})
                % metPubChemID
%                model.metPubChemID(c) = chl.metPubChemID(h);
%            end
%            if isempty(model.metInChIString{c})
                % metInChIString
%                model.metInChIString(c) = chl.metInChIString(h);
%            end
%            if model.metCharge(c) == 0
                % metCharge
%                model.metCharge(c) = chl.metCharge(h);
%            end
        end
        % S
        model.S(c,n) = metMatrix(i,2);
    end
    
    % Handeling genes
%    rxGenes = find(chl.rxnGeneMat(j,:) ~= 0);
%    rxGeneNames = chl.genes(rxGenes);
%    x = length(model.genes) + 1;
%    if ~isempty(rxGenes)
%        for b = 1:length(rxGenes)
%            if isempty(smatch(model.genes, char(rxGeneNames(b)), 'exact')) % Gene is not already found in model
%                % genes
%                model.genes(x) = rxGeneNames(b);
%                % rxnGeneMat
%                model.rxnGeneMat(n,x) = 1;
%                x = x + 1;
%            elseif ~isempty(smatch(model.genes, char(rxGeneNames(b)), 'exact')) % Gene already exists in model
%               if length(smatch(model.genes, char(rxGeneNames(b)), 'exact')) == 1
%                    yPos = smatch(model.genes, char(rxGeneNames(b)), 'exact');
%                    % rxnGeneMat
%                    model.rxnGeneMat(n,yPos) = 1;
%                elseif length(smatch(model.genes, char(rxGeneNames(b)), 'exact')) > 1
%                    warning('Multiple gene hits found for gene %s\n', char(rxGeneNames(b)))
%                    disp(' ')
%                end
%            end
%        end
%    elseif isempty(rxGenes)
%        currentSize = size(model.rxnGeneMat,2);
%        model.rxnGeneMat(n,:) = zeros(1,currentSize);
%    end
end

% Add new description
cytDescription = cyt.description;
chlDescription = chl.description;
model.description = sprintf('Merged model made from %s and %s', cytDescription, chlDescription);

% Getting rid of border reactions

borderRxns = smatch(model.rxns, 'B_');
for j = length(borderRxns):-1:1
    model.confidenceScores(borderRxns(j)) = [];
    model = removeRxns(model, model.rxns(borderRxns(j)));
end

clear i j n metMatrix metIDs m met h c rxGenes rxGeneNames x b yPos currentSize cytDescription chlDescription borderRxns modelMetIndex

%% 6 - Adjusting biomass

% Identifying biomass reaction and biomass metabolites

fprintf('\nWhen the chloroplast is adding to the overall biomass production, certain biomass metabolites will have a higher production rate\n')
adjust = input('Upscale biomass metabolites produced by the chloroplast by 20 percent? (1 = yes, 0 = no): ');
if adjust == 1
    bioRxns = smatch(model.rxns, 'bio');
    rxList = ReactionNames(model);
    if length(bioRxns) == 1
        fprintf('One possible biomass reaction found!\n')
        fprintf('ID:\n   %s\n', char(model.rxns(bioRxns)))
        fprintf('Name:\n   %s\n', char(model.rxnNames(bioRxns)))
        fprintf('Reaction equation:\n   %s\n\n', char(rxList(bioRxns)))
        correct = input('Is this the correct biomass equation? (yes = 1, no = 0): ');
        disp(' ')
        if correct == 1
            bioRxn = bioRxns;
            keepSearching = 0;
            clear bioRxns
        elseif correct == 0
            keepSearching = 1;
        end
        clear correct
    elseif length(bioRxns) > 1
        for i = 1:length(bioRxns)
            fprintf('**** Search result %i ****\n', i)
            fprintf('ID:\n   %s\n', char(model.rxns(bioRxns(i))))
            fprintf('Name:\n   %s\n', char(model.rxnNames(bioRxns(i))))
            fprintf('Reaction equation:\n   %s\n\n', char(rxList(bioRxns(i))))
        end
        indx = input('Choose index of correct biomass reaction (if none of these is the correct biomass reaction, hit 0): ');
        disp(' ')
        if indx > 0
            bioRxn = bioRxns(indx);
            keepSearching = 0;
            clear bioRxns
        else
            keepSearching = 1;
        end    
    elseif isempty(bioRxns)
        keepSearching = 1;    
    end
    while keepSearching == 1
        searchTerm = input('Choose new search term for biomass reaction: ', 's');
        disp(' ')
        list = input('Search in rxn ID list (1) or rxn name list (2)? ');
        disp(' ')
        if list == 1
            res = smatch(model.rxns, searchTerm);
        elseif list == 2
            res = smatch(model.rxnNames, searchTerm);
        end
        if ~isempty(res)
            for i = 1:length(res)
                fprintf('**** Search result %i ****\n', i)
                fprintf('ID:\n   %s\n', char(model.rxns(res(i))))
                fprintf('Name:\n   %s\n', char(model.rxnNames(res(i))))
                fprintf('Reaction equation:\n   %s\n\n', char(rxList(res(i))))
            end
        else
            fprintf('No search results to show\n\n')
        end
        indx = input('Index of correct reaction (if none of these is the correct biomass reaction, hit 0): ');
        disp(' ')
        if indx > 0
            bioRxn = res(indx);
            keepSearching = 0;
        end
    end
    clear keepSearching
    
    % Determining biomass metabolites
    
    bioReactants = find(model.S(:,bioRxn) < 0);
    bioProducts = find(model.S(:,bioRxn) > 0);
    bioReacID = model.mets(bioReactants);
    bioProdID = model.mets(bioProducts);
    
    % Going through biomass reactants to see if chloroplast version exists. If yes, finding out if it is transported accross chloroplast membrane. If it is, and if transport reaction is either irreversible out of chloroplast or reversible, upscaling S-coefficient by 10
    
    for i = 1:length(bioReacID)
        bioMet = char(bioReacID(i));
        hMet = sprintf('%s[h]', bioMet(1:end-3));
        hit = smatch(model.mets, hMet, 'exact');
        if ~isempty(hit) && ~strcmp(char(model.metKEGGID(bioReactants(i))), 'C00001') % Chloroplast metabolite exists and is not water (will be dealt with later in script)
            bioMetRxns = find(model.S(bioReactants(i),:) ~= 0); % Finding reactions biomass metabolite is involved in
            hMetRxns = find(model.S(hit,:) ~= 0); % Finding reaction chloroplast version of biomass metabolite is involved in
            transportRxns = intersect(bioMetRxns, hMetRxns); % Finding transport reactions
            if ~isempty(transportRxns) % Transport reaction(s) exist(s)
                upscale = 0; % Marker for upscaling stoichiometric coefficient in biomass
                for j = 1:length(transportRxns)
                    if (model.rev(transportRxns(j)) == 0) && (model.S(hit,transportRxns(j)) < 0) % Irreversible transport of chloroplast metabolite out of chloroplast
                        upscale = 1;
                    elseif model.rev(transportRxns(j)) == 1
                        upscale = 1;
                    end
                end
                if upscale == 1 % Upscaling S-coefficient by 10 %
                    currentScoeff = nonzeros(model.S(bioReactants(i),bioRxn));
                    newScoeff = currentScoeff * 1.1;
                    metKEGGID = char(model.metKEGGID(bioReactants(i)));
                    if strcmp(metKEGGID, 'C00002') % Saving S-coefficient for ATP
                        oATPcoeff = currentScoeff;
                        nATPcoeff = newScoeff;
                    end
                    if strcmp(metKEGGID, 'C00003') % Saving S-coefficient for NAD+
                        oNADcoeff = currentScoeff;
                        nNADcoeff = newScoeff;
                    end
                    if strcmp(metKEGGID, 'C00005') % Saving S-coefficient for NADPH
                        oNADPHcoeff = currentScoeff;
                        nNADPHcoeff = newScoeff;
                    end
                    if strcmp(metKEGGID, 'C00024') % Saving S-coefficient for Acetyl-CoA
                        oAcCoAcoeff = currentScoeff;
                        nAcCoAcoeff = newScoeff;
                    end
                    model.S(bioReactants(i),bioRxn) = newScoeff;
                end
            end
        elseif strcmp(model.metKEGGID(bioReactants(i)), 'C00001')
            H2Onr = i;
        end
    end
    clear bioMet hMet hit bioMetRxns hMetRxns transportRxns upscale i j currentScoeff newScoeff metKEGGID
    
    if exist('H2Onr', 'var') % Treating water separately
        oH2Ocoeff = nonzeros(model.S(H2Onr,bioRxn));
        nH2Ocoeff = oH2Ocoeff + (nATPcoeff-oATPcoeff);
        if abs(oH2Ocoeff) == abs(oATPcoeff)
            model.S(bioReactants(H2Onr),bioRxn) = nH2Ocoeff;
        else
            fprintf('Stoichiometric coefficient for water is usually tied to the stoichiometric coefficitne for ATP\nbut does not seem to be so in this case\n')
            fprintf('Other biomass metabolites involved with the chloroplast has been upscaled by 10 percent\n')
            upscale = input('Upscale stoichiometric coefficient for water by 10 % (yes = 1, no = 0): ');
            if upscale == 1
                model.S(bioReactants(index),bioRxn) = (oH2Ocoeff * 1.1);
            else
                fprintf('Old stoichiometric coefficient: %f\n', oH2Ocoeff)
                model.S(bioReactants(index),bioRxn) = input('Choose new stoichiometric coefficient for water: ');
                disp(' ')
            end
        end
    end
    clear index H2Onr oH2Ocoeff nH2Ocoeff bioMet hMet hit
    
    % Going through products of biomass reactions
    
    for i = 1:length(bioProdID)
        metKEGGID = char(model.metKEGGID(bioProducts(i)));
        metName = char(model.metNames(bioProducts(i)));
        if strcmp(metKEGGID, 'C00008') % Treating ADP
            if exist('oATPcoeff', 'var')
                oADPcoeff = nonzeros(model.S(bioProducts(i),bioRxn));
                if oADPcoeff == abs(oATPcoeff)
                    nADPcoeff = oADPcoeff + abs(nATPcoeff-oATPcoeff);
                    model.S(bioProducts(i),bioRxn) = nADPcoeff;
                else
                    fprintf('Stoichiometric components of ATP and ADP does not match. Please define manually\n')
                    fprintf('Old stoichiometric coefficient for ATP: %f\n', oATPcoeff)
                    fprintf('New stoichiometric coefficient for ATP: %f\n', nATPcoeff)
                    fprintf('Old stoichiometric coefficient for ADP: %f\n', oADPcoeff)
                    nADPcoeff = input('Define new stoichiometric coefficient for ADP: ');
                    disp(' ')
                    model.S(bioProducts(i),bioRxn) = nADPcoeff;
                end
            else
                fprintf('ADP is present as a biomass component, but not ATP. Please check biomass function manually\n')
                fprintf('Stoichiometric coefficient for ADP will not be changed\n\n')
            end
        elseif strcmp(metKEGGID, 'C00080') % Treating protons
            if exist('oATPcoeff', 'var')
                oHcoeff = nonzeros(model.S(bioProducts(i),bioRxn));
                if oHcoeff == abs(oATPcoeff)
                    nHcoeff = oHcoeff + abs(nATPcoeff-oATPcoeff);
                    model.S(bioProducts(i),bioRxn) = nHcoeff;
                else
                    fprintf('Stoichiometric components of H+ and ATP does not match.\nThese are usually connected. Please define manually\n')
                    fprintf('Old stoichiometric coefficient for ATP: %f\n', oATPcoeff)
                    fprintf('New stoichiometric coefficient for ATP: %f\n', nATPcoeff)
                    fprintf('Old stoichiometric coefficient for H+: %f\n', oHcoeff)
                    nHcoeff = input('Define new stoichiometric coefficient for H+: ');
                    disp(' ')
                    model.S(bioProducts(i),bioRxn) = nHcoeff;
                end
            else
                fprintf('H+ is present as a biomass component, but not ATP. These are usually connected.\nPlease check biomass function manually\n')
                fprintf('Stoichiometric coefficient for H+ will not be changed\n\n')
            end
        elseif strcmp(metKEGGID, 'C00009') % Treating inorganic phosphate Pi
            if exist('oATPcoeff', 'var')
                oPIcoeff = nonzeros(model.S(bioProducts(i),bioRxn));
                if oPIcoeff == abs(oATPcoeff)
                    nPIcoeff = oPIcoeff + abs(nATPcoeff-oATPcoeff);
                    model.S(bioProducts(i),bioRxn) = nPIcoeff;
                else
                    fprintf('Stoichiometric components of inorganic phosphate (Pi) and ATP does not match.\nThese are usually connected. Please define manually\n')
                    fprintf('Old stoichiometric coefficient for ATP: %f\n', oATPcoeff)
                    fprintf('New stoichiometric coefficient for ATP: %f\n', nATPcoeff)
                    fprintf('Old stoichiometric coefficient for Pi: %f\n', oPIcoeff)
                    nPIcoeff = input('Define new stoichiometric coefficient for Pi: ');
                    disp(' ')
                    model.S(bioProducts(i),bioRxn) = nPIcoeff;
                end
            else
                fprintf('Inorganic phosphate (Pi) is present as a biomass component, but not ATP.\nThese are usually connected. Please check biomass function manually\n')
                fprintf('Stoichiometric coefficient for Pi will not be changed\n\n')
            end
        elseif strcmp(metKEGGID, 'C00010') % Treating Coenzyme A
            if exist('oAcCoAcoeff', 'var')
                oCOAcoeff = nonzeros(model.S(bioProducts(i),bioRxn));
                if oCOAcoeff == abs(oAcCoAcoeff)
                    nCOAcoeff = abs(nAcCoAcoeff);
                    model.S(bioProducts(i),bioRxn) = nCOAcoeff;
                else
                    fprintf('Stoichiometric components of Coenzyme A and Acetyl-CoA does not match.\nThese are usually connected. Please define manually\n')
                    fprintf('Old stoichiometric coefficient for Acetyl-CoA: %f\n', oAcCoAcoeff)
                    fprintf('New stoichiometric coefficient for Acetyl-CoA: %f\n', nAcCoAcoeff)
                    fprintf('Old stoichiometric coefficient for Coenzyme A: %f\n', oCOAcoeff)
                    nCOAcoeff = input('Define new stoichiometric coefficient for Coenzyme A: ');
                    disp(' ')
                    model.S(bioProducts(i),bioRxn) = nCOAcoeff;
                end
            else
                fprintf('Coenzyme A is present as a biomass component, but not Acetyl-CoA.\nThese are usually connected. Please check biomass function manually\n')
                fprintf('Stoichiometric coefficient for Coenzyme A will not be changed\n\n')
            end
        elseif strcmp(metKEGGID, 'C00004') % NADH
            if exist('oNADcoeff', 'var')
                oNADHcoeff = nonzeros(model.S(bioProducts(i),bioRxn));
                if oNADHcoeff == abs(oNADcoeff)
                    nNADHcoeff = abs(nNADcoeff);
                    model.S(bioProducts(i),bioRxn) = nNADHcoeff;
                else
                    fprintf('Stoichiometric components of NADH and NAD+ does not match.\nThese are usually connected. Please define manually\n')
                    fprintf('Old stoichiometric coefficient for NAD+: %f\n', oNADcoeff)
                    fprintf('New stoichiometric coefficient for NAD+: %f\n', nNADcoeff)
                    fprintf('Old stoichiometric coefficient for NADH: %f\n', oNADHcoeff)
                    nNADHcoeff = input('Define new stoichiometric coefficient for NADH: ');
                    disp(' ')
                    model.S(bioProducts(i),bioRxn) = nNADHcoeff;
                end
            elseif isempty(smatch(model.metKEGGID(bioReactants), 'C00003'))
                fprintf('NADH is present as a biomass component, but not NAD+.\nThese are usually connected. Please check biomass function manually\n')
                fprintf('Stoichiometric coefficient for NADH will not be changed\n\n')
            end
        elseif strcmp(metKEGGID, 'C00006') % NADP+
            if exist('oNADPHcoeff', 'var')
                oNADPcoeff = nonzeros(model.S(bioProducts(i),bioRxn));
                if oNADPcoeff == abs(oNADPHcoeff)
                    nNADPcoeff = abs(nNADPHcoeff);
                    model.S(bioProducts(i),bioRxn) = nNADPcoeff;
                else
                    fprintf('Stoichiometric components of NADP+ and NADPH does not match.\nThese are usually connected. Please define manually\n')
                    fprintf('Old stoichiometric coefficient for NADPH: %f\n', oNADPHcoeff)
                    fprintf('New stoichiometric coefficient for NADPH: %f\n', nNADPHcoeff)
                    fprintf('Old stoichiometric coefficient for NADP+: %f\n', oNADPcoeff)
                    nNADPcoeff = input('Define new stoichiometric coefficient for NADP+: ');
                    disp(' ')
                    model.S(bioProducts(i),bioRxn) = nNADPcoeff;
                end
            elseif isempty(smatch(model.metKEGGID(bioReactants), 'C00005'))
                fprintf('NADP+ is present as a biomass component, but not NADPH.\nThese are usually connected. Please check biomass function manually\n')
                fprintf('Stoichiometric coefficient for NADP+ will not be changed\n\n')
            end
        else
            fprintf('Manual upscaling of metabolite %s:\n\n', metName(1:end-4))
            fprintf('Upscaled stoichiometric coefficients for the reactant side of the biomass function:\n')
            disp([model.metNames(bioReactants), num2cell(nonzeros(model.S(bioReactants,bioRxn)))])
            disp(' ')
            fprintf('Current stoichiometric coefficient for metabolite %s: %f\n\n', metName(1:end-4), nonzeros(model.S(bioProducts(i),bioRxn)))
            replace = input('Enter new stoichiometric coefficient or keep the old one? (enter new = 1, keep old = 0): ');
            disp(' ')
            if replace == 1
                nScoeff = input('New stoichiometric coefficient: ');
                disp(' ')
                model.S(bioProducts(i),bioRxn) = nScoeff;
            end
        end
    end
    clear metKEGGID H2Onr nAcCoAcoeff oAcCoAcoeff nADPcoeff oADPcoeff nATPcoeff oATPcoeff nCOAcoeff oCOAcoeff nH2Ocoeff oH2Ocoeff nHcoeff oHcoeff nNADcoeff oNADcoeff nNADHcoeff oNADHcoeff nNADPcoeff oNADPcoeff nNADPHcoeff oNADPHcoeff nPIcoeff oPIcoeff replace nScoeff
    
    % Adding biomass components from chloroplast
    
    chlBioMets = smatch(model.metNames, 'Chloroplast'); % Gathering list of metabolites that should be added to biomass function
    for i = length(chlBioMets):-1:1
        metName = char(model.metNames(chlBioMets(i)));
        if ~strcmp(metName(1:11), 'Chloroplast')
            chlBioMets(i) = [];
        end
    end
    for i = length(chlBioMets):-1:1
        metName = char(model.metNames(chlBioMets(i)));
        if isempty(strfind(metName, '[c]'))
            chlBioMets(i) = [];
        end
    end
    clear i metName
    
    bioMets = [bioReactants; bioProducts];
    clear bioProdID bioProducts bioReacID bioReactants
    for i = 1:length(chlBioMets)
        fprintf('* Adding metabolite %s to biomass function *\n\n', char(model.metNames(chlBioMets(i))))
        dispHelp = input('Display other stoichiometric coefficients in biomass for help? (yes = 1, no = 0): ');
        disp(' ')
        if dispHelp == 1
            disp([model.metNames(bioMets), num2cell(nonzeros(model.S(bioMets,bioRxn)))])
        end
        sCoeff = input('Choose stoichiometric coefficient for metabolite: ');
        disp(' ')
        disp(' ')
        if sCoeff > 0
            sCoeff = sCoeff * (-1);
        end
        model.S(chlBioMets(i),bioRxn) = sCoeff;
    end
    clear bioMets dispHelp sCoeff
    clc
end

clc
fprintf('Congratulations, you successfully merged the cytosol model and the chloroplast model!\n\n')
load handel.mat
sound(y,Fs)
clear y Fs