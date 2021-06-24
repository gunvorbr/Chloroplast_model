function finalResult = findMet(model, met)

% Searches model for metabolites with a certain name, or where the
% metabolite name includes a certain phrase (ex. 'phenyl')
% Displays number, met-ID, metName and KEGG-ID for search hits

%   by Gunvor Røkke, NTNU, 2021

result = cell(1,4); % Make space for #, ID, name and KEGG
n = 1;

for i = 1:length(model.metNames) % Searches metNames vector
    if strfind(lower(char(model.metNames(i))), lower(char(met)))
        result(n,1) = num2cell(i);
        result(n,2) = model.mets(i);
        result(n,3) = model.metNames(i);
        result(n,4) = model.metKEGGID(i);
        n = n + 1;
    end
end

for i = 1:length(model.mets) % Searches mets vector
    if strfind(lower(char(model.mets(i))), lower(char(met)))
        result(n,1) = num2cell(i);
        result(n,2) = model.mets(i);
        result(n,3) = model.metNames(i);
        result(n,4) = model.metKEGGID(i);
        n = n + 1;
    end
end

for i = 1:length(model.mets) % Searches metKEGGID vector
    if strfind(lower(char(model.metKEGGID(i))), lower(char(met)))
        result(n,1) = num2cell(i);
        result(n,2) = model.mets(i);
        result(n,3) = model.metNames(i);
        result(n,4) = model.metKEGGID(i);
        n = n + 1;
    end
end

finalResult = cell(1,4);
k = 1;

for i = 1:length(model.mets)
    if find(cell2mat(result(:,1)) == i)
        hits = find(cell2mat(result(:,1)) == i);
        finalResult(k,1) = result(hits(1),1);
        finalResult(k,2) = result(hits(1),2);
        finalResult(k,3) = result(hits(1),3);
        finalResult(k,4) = result(hits(1),4);
        k = k + 1;
    end
end