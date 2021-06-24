function rxns = KEGGrxGet(model)

% Input:
% model -> Is being used to search for metabolites
% KEGG-ids of metabolites are being used to search KEGG for common
% reactions between metabolites

%   by Gunvor Røkke, NTNU, 2021

go = 1;
KEGGlist = cell(1,1);
m = 1;

while go == 1
    met = input('Metabolite? ', 's');
    disp(' ')
    res = findPattern(model.metNames, met);
    if cellfun(@isempty, res(1,:)) == [1,1] % Results for metabolite search not found
        disp('Metabolite not found!')
        disp(' ')
        KEGGlist(m,1) = cellstr(input('Check KEGG ID for metabolite manually: ', 's'));
    elseif cellfun(@isempty, res(1,:)) ~= [1,1] % Results for metabolite search found
        disp([transpose(num2cell(1:size(res,1))), res])
        pick = input('Correct metabolite # (if not in list, choose 0): ');
        if pick == 0 % Metabolite not found
            KEGGlist(m,1) = cellstr(input('KEGG ID for metabolite: ', 's'));
        elseif pick ~= 0 % Metabolite found in model
            nr = cell2mat(res(pick,1));
            KEGGlist(m,1) = model.metKEGGID(nr);
        end
    end
    m = m + 1;
    go = input('Continue searching for metabolites? (1 = yes, 0 = no): ');
    disp(' ')
end

url = 'http://rest.kegg.jp/get/';

for i = 1:length(KEGGlist)
    res = urlread(strcat(url, char(KEGGlist(i))));
    if ~isempty(res)
        start = strfind(res, 'REACTION') + 12;
        if strfind(res, 'PATHWAY')
            slutt = strfind(res, 'PATHWAY') - 2;
        elseif strfind(res, 'MODULE')
            slutt = strfind(res, 'MODULE') - 2;
        elseif strfind(res, 'ENZYME')
            slutt = strfind(res, 'ENZYME') - 2;
        elseif strfind(res, 'BRITE')
            slutt = strfind(res, 'BRITE') - 2;
        elseif strfind(res, 'DBLINKS')
            slutt = strfind(res, 'DBLINKS') - 2;
        end
        if ~exist('rxns')
            rxns = strsplit(res(start:slutt));
        elseif exist('rxns')
            newRxns = strsplit(res(start:slutt));
            rxns = intersect(newRxns, rxns);
        end
    elseif isempty(res)
        disp('Uh oh! Something went terribly wrong :-(')
    end
end

rxns = transpose(rxns);
url = 'http://rest.kegg.jp/get/';

disp(' ')
disp(' ')
disp(' ')

for i = 1:length(rxns)
    fprintf('*** %i: Reaction %s ***', i, char(rxns(i)))
    disp(' ')
    res = urlread(strcat(url, char(rxns(i))));
    start = strfind(res, 'DEFINITION') + 12;
    slutt = strfind(res, 'EQUATION') - 2;
    disp(res(start:slutt))
    disp(' ')
    disp(' ')
    disp(' ')
end

disp('The end :-)')