function model = fillKEGGIDholes(model)

% Filling holes in KEGG vector that was not filled by metKEGGIDsearch.
% Allows for modification of search term
%   Gunvor Røkke, NTNU, 2021

adresse = 'http://rest.kegg.jp/find/compound/';

disp(' ')
for i = 1:length(model.metKEGGID)
    if isempty(model.metKEGGID{i})
        search = 1;
        while search == 1
            fprintf('*** Metabolite %i  of %i: %s (%s) (%s) ***\n\n', i, length(model.metNames), char(model.metNames(i)), char(model.mets(i)), char(model.metFormulas(i)))
            metSok = input('Search term: ', 's');
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
            
            disp(' ')
            if strcmp(char(tabell(1,1)), 'NaN')
                disp('No hits!')
            else
                disp(tabell)
            end
            
            nr = input('Index of correct metabolite (metabolite not found, or new search required -> press 0): ');
            disp(' ')
            
            if nr > 0
                search = 0;
            else
                again = input('Search again? (yes = 1, no = 0): ');
                disp(' ')
                if again == 0
                    search = 0;
                end
            end
        end
        
        if nr > 0
            empty = cellfun(@isempty, tabell(:,2));
            for j = 1:length(empty)
                if empty(j) == 1
                    tabell(j,2) = cellstr('0');
                elseif empty(j) == 0
                    tabell(j,2) = cellstr(num2str(cell2mat(tabell(j,2))));
                end
            end
            row = strmatch(num2str(nr), tabell(:,2), 'exact');
            model.metKEGGID(i) = tabell(row,1);
        end
        dispMet(model,i)
        clear empty j nr row tabell
        disp(' ')
        disp(' ')
    end
end
clc