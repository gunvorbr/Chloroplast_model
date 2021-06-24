function chl = changeOrganismMode(chl)

% Changes organism mode in model with multiple organism modes
%   by Gunvor Røkke, NTNU, 2021

disp(' ')
org = input('Code letter for organism? ', 's');

totalOrgRxns = double.empty(0,1);
n = 1;
for i = 1:length(chl.rxns)
    if strfind(char(chl.rxns(i)), '@')
        totalOrgRxns(n,1) = i;
        n = n + 1;
    end
end

relevantOrgRxns = double.empty(0,1);
n = 1;
for a = 1:length(totalOrgRxns)
    i = totalOrgRxns(a);
    if strfind(char(chl.rxns(i)), sprintf('@%s', org))
        relevantOrgRxns(n,1) = i;
        n = n + 1;
    end
end

for a = 1:length(totalOrgRxns)
    i = totalOrgRxns(a);
    if isempty(find(relevantOrgRxns == i)) % Reaction belongs to different organism
        chl.ub(i) = 0;
        chl.lb(i) = 0;
    elseif ~isempty(find(relevantOrgRxns == i)) % Reaction belongs to relevant organism
        chl.ub(i) = 1000;
        if chl.rev(i) == 1
            chl.lb(i) = -1000;
        end
    end
end

% Changing reaction to be optimized

if isempty(strfind(chl.description, 'Merged'))
    makeZero = find(chl.c ~= 0);
    for i = 1:length(makeZero)
        chl.c(makeZero(i)) = 0;
    end
    
    optimRxn = sprintf('%s_bio', org);
    bioRxn = smatch(chl.rxns, optimRxn);
    
    chl.c(bioRxn) = 1;
end
        