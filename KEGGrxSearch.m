function rxns = KEGGrxSearch(KEGGlist)

% Accepts a list of metabolite KEGG-IDs as input
% Searches KEGG for reactions the given metabolit(s) participate in
% Prints rx equations and EC-numbers
% Outputs list of rx KEGG-IDs

%   by Gunvor Røkke, NTNU, 2021

if ischar(KEGGlist)
    KEGGlist = cellstr(KEGGlist);
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
        warning('Uh oh! Something went terribly wrong :-(')
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
    disp(' ')
    res = urlread(strcat(url, char(rxns(i))));
    start = strfind(res, 'DEFINITION') + 12;
    slutt = strfind(res, 'EQUATION') - 2;
    disp(res(start:slutt))
    disp(' ')
    
    eqStart = strfind(res, 'EQUATION') + 12;
    eqSlutt = eqStart;
    while isempty(strfind(res(eqSlutt), sprintf('\n')))
        eqSlutt = eqSlutt + 1;
    end
    disp(res(eqStart:eqSlutt-1))
    disp(' ')
    if strfind(res, 'ENZYME')
        eStart = strfind(res, 'ENZYME') + 12;
        if strfind(res, 'PATHWAY')
            eSlutt = strfind(res, 'PATHWAY') - 2;
        elseif strfind(res, 'MODULE')
            eSlutt = strfind(res, 'MODULE') - 2;
        elseif strfind(res, 'ORTHOLOGY')
            eSlutt = strfind(res, 'ORTHOLOGY') - 2;
        elseif strfind(res, '///')
            eSlutt = strfind(res, '///') - 2;
        end
        disp(res(eStart:eSlutt))
    end
    disp(' ')
    disp(' ')
    disp(' ')
end

disp('The end :-)')