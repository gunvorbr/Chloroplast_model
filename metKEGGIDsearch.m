function model = metKEGGIDsearch(model)

% Runs through metabolite names in model, and tries to find corresponding
% metKEGGID
% Lets you choose correct ID
% Stores metabolite KEGG-IDs in model.metKEGGID

%   by Gunvor Røkke, NTNU, 2021

adresse = 'http://rest.kegg.jp/find/compound/';

% Creating array of root metabolite IDs
rootMets = cell(length(model.mets),1);
metTest = char(model.mets(1));
if strcmp(metTest(end-2), '[')
    backSteps = 3;
elseif strcmp(metTest(end-1), '_')
    backSteps = 2;
elseif strcmp(metTest(end-2), '_')
    backSteps = 3;
end
for i = 1:length(model.mets)
    met = char(model.mets(i));
    rootMets(i) = cellstr(met(1:end-backSteps));
end
clc

for i = 1:length(model.metNames)
    fprintf('*** Metabolite %i of %i: %s ***\n\n', i, length(model.mets), char(model.metNames(i)))
    if isempty(model.metKEGGID{i})
        % Checking if metabolite has already been searched for
        parallelMets = smatch(rootMets, char(rootMets(i)), 'exact');
        parallelKEGGs = model.metKEGGID(parallelMets);
        KEGGID = 'X';
        for k = 1:length(parallelKEGGs)
            if ~isempty(parallelKEGGs{k})
                KEGGID = char(parallelKEGGs(k));
            end
        end
        if ~strcmp(KEGGID, 'X')
            model.metKEGGID(i) = cellstr(KEGGID);
        else
            metSok = char(model.metNames(i));
            if length(metSok) > 5
                if strcmp(metSok(end-3:end-2), ' [')
                    metSok = metSok(1:end-4);
                end
            end
            metSok = strrep(metSok, ' ', '%20');
            metSok = strrep(metSok, ':', '%3A');
            metSok = strrep(metSok, ';', '%3B');
            metSok = strrep(metSok, ',', '%2C');
            metSok = strrep(metSok, '(', '%28');
            metSok = strrep(metSok, ')', '%29');
            metSok = strrep(metSok, '[', '%5B');
            metSok = strrep(metSok, ']', '%5D');
            metSok = strrep(metSok, '-', '%20');
            metSok = strrep(metSok, '/', '%20');
            
            res = urlread(strcat(adresse, metSok));
            
            if length(res) == 1
                tabell = cellstr('NaN');
            else
                rows = strfind(res, sprintf('\n'));
                mTabell = cell(length(rows),1);
                mTabell(1,1) = cellstr(res(5:rows(1)-1));
                if length(rows) > 1
                    for j = 1:(length(rows)-1)
                        mTabell(j+1,1) = cellstr(res(rows(j)+5:rows(j+1)-1));
                    end
                end
                clear j
                tabell = cell(1,3);
                t = 1;
                for j = 1:length(mTabell)
                    row = char(mTabell(j));
                    semiCols = strfind(row, ';');
                    tabell(t,1) = cellstr(row(1:6));
                    tabell(t,2) = num2cell(j);
                    if isempty(semiCols)
                        tabell(t,3) = cellstr(row(8:end));
                        t = t + 1;
                    elseif length(semiCols) == 1
                        tabell(t,3) = cellstr(row(8:semiCols-1));
                        t = t + 1;
                        tabell(t,3) = cellstr(row(semiCols+2:end));
                        t = t + 1;
                    else
                        tabell(t,3) = cellstr(row(8:semiCols(1)-1));
                        t = t + 1;
                        for k = 1:(length(semiCols)-1)
                            tabell(t,3) = cellstr(row(semiCols(k)+2:semiCols(k+1)-1));
                            t = t + 1;
                        end
                        tabell(t,3) = cellstr(row(semiCols(end)+2:end));
                        t = t + 1;
                    end
                end
            end
            clear j metSok mTabell res row rows semiCols t
            
            if ~strcmp(char(tabell(1,1)), 'NaN')
                if length(cell2mat(tabell(:,2))) == 1 % Choosing KEGG-ID automatically if only one hit
                    model.metKEGGID(i) = tabell(1,1);
                else
                    % Checking if exact hit is found in table
                    metName = char(model.metNames(i));
                    metName = metName(1:end-4);
                    metName = strrep(metName, '-', ' ');
                    tableHits = double.empty(0,1);
                    n = 1;
                    for k = 1:size(tabell,1)
                        tableStr = char(tabell(k,3));
                        tableStr = strrep(tableStr, '-', ' ');
                        if strcmpi(metName, tableStr)
                            tableHits(n,1) = k;
                            n = n + 1;
                        end
                    end
                    if length(tableHits) == 1
                        k = tableHits;
                        while isempty(tabell{k,1})
                            k = k - 1;
                        end
                        model.metKEGGID(i) = tabell(k,1);
                    else
                        fprintf('*** Metabolite %i out of %i, %s (%s) ***\n', i, length(model.mets), char(model.metNames(i)), char(model.mets(i)))
                        disp(' ')
                        if strcmp(char(tabell(1,1)), 'NaN')
                            disp('Metabolite not found in KEGG')
                        else
                            disp(tabell)
                            
                            nr = input('Number of correct KEGG ID (if not in list, hit 0): ');
                            if nr ~= 0
                                empty = cellfun(@isempty, tabell(:,2));
                                for j = 1:length(empty)
                                    if empty(j) == 1
                                        tabell(j,2) = cellstr('0');
                                    elseif empty(j) == 0
                                        tabell(j,2) = cellstr(num2str(cell2mat(tabell(j,2))));
                                    end
                                end
                                
                                row = smatch(tabell(:,2), num2str(nr), 'exact');
                                model.metKEGGID(i) = tabell(row,1);
                            end
                        end
                    end
                end
            end
            clear empty j nr row tabell
            disp(' ')
        end
        if ~isempty(model.metKEGGID{i})
            fprintf('KEGGID set to %s\n', char(model.metKEGGID(i)))
        end
        disp(' ')
        disp(' ')
    end
end
clc