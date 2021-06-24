function cyt = KEGGIDfirstAid(cyt, chl)

% Attempts to find KEGG-ID for all metabolites in a model lacking KEGG-IDs
%   by Gunvor Røkke, NTNU, 2021

chlcMets = double.empty(0,1);
n = 1;

disp(' ')
for i = 1:size(chl.mets) % Looping through chl metabolites, if metabolite is in c-compartment, adding to chlcMets
    currentMet = char(chl.mets(i));
    if ~isempty(strfind(currentMet, '[c]'))
        chlcMets(n,1) = i;
        n = n + 1;
    end
end
clear currentMet i n

chlcMets = unique(chlcMets);
chlcNames = chl.metNames(chlcMets);

for i = 1:length(chlcNames)
    fprintf('(Met %i of %i)\n\n', i, length(chlcMets))
    currentName = char(chlcNames(i));
    searchTerm = currentName(1:end-4);
    if isempty(strfind(lower(searchTerm), 'chloroplast'))
        hits = findPattern(cyt.metNames, searchTerm); % Searching cyt model for current chlcmet
        if ~isempty(hits{1,1}) % Hits exist in cyt for current metabolite
            exactHits = double.empty(0,1);
            n = 1;
            for k = 1:size(hits,1) % Checking if exact hits exist
                currentHit = char(hits(k,2));
                if strcmpi(searchTerm, currentHit(1:end-4))
                    exactHits(n,1) = cell2mat(hits(k,1));
                    n = n + 1;
                end
            end
            if ~isempty(exactHits)
                for k = 1:length(exactHits)
                    if strfind(char(cyt.metNames(exactHits(k))), '[c]')
                        exactHit = exactHits(k);
                    end
                end
                if ~exist('exactHit', 'var')
                    noExactHits = 1;
                else
                    noExactHits = 0;
                end
            elseif isempty(exactHits)
                noExactHits = 1;
            end
            if noExactHits == 1
                fprintf('Hits in cytosol model for metabolite %s (%s):\n\n', searchTerm, char(chl.metKEGGID(chlcMets(i))))
                disp([num2cell(transpose(1:size(hits,1))), hits])
                disp(' ')
                numb = input('Index of correct metabolite? (press 0 if correct metabolite is not present): ');
                disp(' ')
                if numb > 0
                    exactHit = cell2mat(hits(numb,1));
                elseif numb == 0
                    again = input('Search again with modified search term? (yes = 1, no = 0): ');
                    disp(' ')
                    while again == 1
                        searchTerm = input('Modified search term: ', 's');
                        disp(' ')
                        hits = findPattern(cyt.metNames, searchTerm);
                        if ~isempty(hits)
                            disp([num2cell(transpose(1:size(hits,1))), hits])
                            disp(' ')
                            numb = input('Number of correct metabolite? (press 0 if correct metabolite is not present): ');
                            disp(' ')
                            if numb == 0
                                again = input('Search again? (yes = 1, no = 0): ');
                                disp(' ')
                            else
                                exactHit = cell2mat(hits(numb,1));
                                again = 0;
                            end
                        elseif isempty(hits)
                            again = input('Search again? (yes = 1, no = 0): ');
                            disp(' ')
                        end
                    end
                end
            end
        else
            fprintf('No hits found for metabolite *** %s (%s) ***\nPlease search again with modified search term\n\n', searchTerm, char(chl.metKEGGID(chlcMets(i))))
            again = 1;
            while again == 1
                searchTerm = input('Modified search term: ', 's');
                disp(' ')
                hits = findPattern(cyt.metNames, searchTerm);
                if ~isempty(hits)
                    disp([num2cell(transpose(1:size(hits,1))), hits])
                    disp(' ')
                    numb = input('Number of correct metabolite? (press 0 if correct metabolite is not present): ');
                    disp(' ')
                    if numb == 0
                        again = input('Search again? (yes = 1, no = 0): ');
                        disp(' ')
                    else
                        exactHit = cell2mat(hits(numb,1));
                        again = 0;
                    end
                elseif isempty(hits)
                    again = input('Search again? (yes = 1, no = 0): ');
                    disp(' ')
                end
            end
        end
        if exist('exactHit', 'var')
            if exactHit > 0
                cyt.metKEGGID(exactHit) = chl.metKEGGID(chlcMets(i));
                fprintf('%s -> %s\n\n', char(cyt.metNames(exactHit)), char(cyt.metKEGGID(exactHit)))
            elseif exactHit == 0
                fprintf('No KEGGID added for metabolite %s\n\n', currentName)
            end
        else
            fprintf('No KEGGID found for metabolite %s\n\n', currentName)
        end
        clear exactHit currentName searchTerm hits exactHits again
    else
        fprintf('Met %i does not need a corresponding KEGG ID in the cytosol model\n\n', i)
    end
    fprintf('---------------\n\n')
end
