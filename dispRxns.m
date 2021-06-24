function rxns = dispRxns(model)

% Identifies and displays all reactions in a model that a certain metabolite is involved in
% Search with KEGG-ID or metabolite name
%   by Gunvor Røkke, NTNU, 2021

disp(' ')
found = 0;
while found == 0
    searchTerm = input('Search for metabolite: ', 's');
    disp(' ')
    if length(searchTerm) == 6
        if ~isnan(str2double(searchTerm(2:end)))
            res = findPattern(model.metKEGGID, searchTerm);
            KEGGID = 1;
        else
            KEGGID = 0;
            searchList = input('Search in ID (i) or name (n) list? ', 's');
            disp(' ')
            if strcmp(searchList, 'i')
                res = findPattern(model.mets, searchTerm);
            elseif strcmp(searchList, 'n')
                res = findPattern(model.metNames, searchTerm);
            end
        end
    else
        KEGGID = 0;
        searchList = input('Search in ID (i) or name (n) list? ', 's');
        disp(' ')
        if strcmp(searchList, 'i')
            res = findPattern(model.mets, searchTerm);
        elseif strcmp(searchList, 'n')
            res = findPattern(model.metNames, searchTerm);
        end
    end
    if isempty(res{1,1}) % No search results
        if KEGGID == 0
            again = input('No results found. Search again? (yes = 1, no = 0): ');
            disp(' ')
            if again == 0
                return
            end
        elseif KEGGID == 1
            fprintf('No search results found\n\n')
            return
        end
    else % Search results found
        disp([transpose(num2cell(1:size(res,1))), res, model.mets(cell2mat(res(:,1)))])
        pick = input('Correct metabolite # (if not in list, choose 0): ');
        if pick == 0
            again = input('Search again? (yes = 1, no = 0)');
            disp(' ')
            if again == 0
                return
            end
        else % Correct metabolite found
            nr = cell2mat(res(pick,1));
            found = 1;
        end
    end
end

rxns = transpose(find(model.S(nr,:) ~= 0));

rxList = ReactionNames(model);

fprintf('\n\n**** Reactions: ****\n\n')

disp(char(rxList(rxns)))
