function tabell = KmetSearch(metSok)

% Input: Name of metabolite (or fraction of name)
% Script searches KEGG for metabolites matching input
% Outputs table with alternative names and KEGG-IDs of metabolites matching the search word

%   by Gunvor Røkke, NTNU, 2021

adresse = 'http://rest.kegg.jp/find/compound/';

metSok = strrep(metSok, ' ', '%20');
metSok = strrep(metSok, ':', '%3A');
metSok = strrep(metSok, ',', '%2C');
metSok = strrep(metSok, '(', '%28');
metSok = strrep(metSok, ')', '%29');
metSok = strrep(metSok, '[', '%5B');
metSok = strrep(metSok, ']', '%5D');

res = urlread(strcat(adresse, metSok));

if length(res) == 1
    disp(' ')
    disp('No results found!')
    disp(' ')
else
    rows = strfind(res, sprintf('\n'));
    
    mTabell = cell(length(rows),1);
    
    mTabell(1,1) = cellstr(res(5:rows(1)-1));
    
    if length(rows) > 1
        for i = 1:(length(rows)-1)
            mTabell(i+1,1) = cellstr(res(rows(i)+5:rows(i+1)-1));
        end
    end
    
    tabell = cell(1,2);
    t = 1;
    
    for i = 1:length(mTabell)
        row = char(mTabell(i));
        semiCols = strfind(row, ';');
        tabell(t,1) = cellstr(row(1:6));
        if isempty(semiCols)
            tabell(t,2) = cellstr(row(8:end));
            t = t + 1;
        elseif length(semiCols) == 1
            tabell(t,2) = cellstr(row(8:semiCols-1));
            t = t + 1;
            tabell(t,2) = cellstr(row(semiCols+2:end));
            t = t + 1;
        else
            tabell(t,2) = cellstr(row(8:semiCols(1)-1));
            t = t + 1;
            for j = 1:(length(semiCols)-1)
                tabell(t,2) = cellstr(row(semiCols(j)+2:semiCols(j+1)-1));
                t = t + 1;
            end
            tabell(t,2) = cellstr(row(semiCols(end)+2:end));
            t = t + 1;
        end
    end
    
    disp(' ')
    disp(tabell)
end