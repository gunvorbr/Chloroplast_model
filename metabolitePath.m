function result = metabolitePath(model, FBAresult)

% Traces a metabolite path
% Lets user choose if the metabolite should be a substrate or a product
%   by Gunvor Røkke, NTNU, 2021

disp(' ')
searchMet = 1;
while searchMet == 1
    metSearch = input('Search for initial metabolite: ', 's');
    disp(' ')
    res = findPattern(model.metNames, metSearch);
    disp([num2cell(transpose(1:size(res,1))), res])
    indx = input('Number of metabolite (or 0, if you want to search again): ');
    disp(' ')
    metNr = cell2mat(res(indx,1));
    if metNr > 0
        searchMet = 0;
    elseif metNr == 0
        searchMet = input('Search again? (yes = 1, no = 0): ');
        if searchMet == 0
            return
        end
    end
end

subSys = input('Display information about subsystems? (yes = 1, no = 0): ');
disp(' ')

% Displays metabolite path

result = cell(1,1);
n = 1;
rxList = ReactionNames(model);
disp(' ')
dir = input('Study metabolite as substrate (s) or product (p)? ', 's');
fprintf('\n-----\n\n')

go = 1;
while go == 1
    rxns = find(model.S(metNr,:) ~= 0);
    if length(rxns) == 1
        if exist('lastRx', 'var')
            if rxns == lastRx
                disp('--- End of path ---')
                disp(' ')
                if strcmp(dir, 'p')
                    result(n,1) = cellstr(sprintf('-| %s', char(model.metNames(metNr))));
                elseif strcmp(dir, 's')
                    result(n,1) = cellstr(sprintf('%s -|', char(model.metNames(metNr))));
                end
                break
            end
        end
    end
    if strcmp(dir, 's')
        fprintf('Metabolite %s is used in these reactions:\n\n', char(model.metNames(metNr)))
    elseif strcmp(dir, 'p')
        fprintf('Metabolite %s is produced in these reactions:\n\n', char(model.metNames(metNr)))
    end
    resultMets = zeros(1,2);
    k = 1;
    for j = 1:length(rxns)
        if strcmp(dir, 'p')
            if (model.S(metNr,rxns(j)) > 0) && (model.rev(rxns(j)) == 0) % Metabolite is product of irreversible reaction
                fprintf('Rx %i\n', rxns(j))
                disp(char(rxList(rxns(j))))
                if nargin == 2
                    fprintf('Flux = %8.3f\n', FBAresult.x(rxns(j)))
                end
                if subSys == 1
                    fprintf('Sub system: %s\n', char(model.subSystems(rxns(j))))
                end
                disp(' ')
                otherMets = find(model.S(:,rxns(j)) < 0);
                for i = 1:length(otherMets)
                    resultMets(k,1) = otherMets(i);
                    resultMets(k,2) = rxns(j);
                    k = k + 1;
                end
            end
        elseif strcmp(dir, 's')
            if (model.S(metNr,rxns(j)) < 0) && (model.rev(rxns(j)) == 0) % Metabolite is substrate of irreversible reaction
                fprintf('Rx %i\n', rxns(j))
                disp(char(rxList(rxns(j))))
                if nargin == 2
                    fprintf('Flux = %8.3f\n', FBAresult.x(rxns(j)))
                end
                if subSys == 1
                    fprintf('Sub system: %s\n', char(model.subSystems(rxns(j))))
                end
                disp(' ')
                otherMets = find(model.S(:,rxns(j)) > 0);
                for i = 1:length(otherMets)
                    resultMets(k,1) = otherMets(i);
                    resultMets(k,2) = rxns(j);
                    k = k + 1;
                end
            end
        end
        if model.rev(rxns(j)) == 1 % Metabolite is involved in reversible reactions
            fprintf('Rx %i\n', rxns(j))
            disp(char(rxList(rxns(j))))
            if nargin == 2
                fprintf('Flux = %8.3f\n', FBAresult.x(rxns(j)))
            end
            if subSys == 1
                fprintf('Sub system: %s\n', char(model.subSystems(rxns(j))))
            end
            disp(' ')
            if model.S(metNr,rxns(j)) > 0
                otherMets = find(model.S(:,rxns(j)) < 0);
            elseif model.S(metNr,rxns(j)) < 0
                otherMets = find(model.S(:,rxns(j)) > 0);
            end
            for i = 1:length(otherMets)
                resultMets(k,1) = otherMets(i);
                resultMets(k,2) = rxns(j);
                k = k + 1;
            end
        end
    end
    if exist('resultMets', 'var')
        if strcmp(dir, 's')
            fprintf('\nMetabolite could turn into these metabolites:\n\n') % Displaying new possible metabolites to follow
        elseif strcmp(dir, 'p')
            fprintf('\nMetabolite could be formed from these metabolites:\n\n')
        end
        for i = 1:size(resultMets,1)
            fprintf('[%i]  %s\n', i, char(model.metNames(resultMets(i,1))))
            fprintf('       in reaction %i, %s\n\n', resultMets(i,2), char(model.rxnNames(resultMets(i,2))))
        end
        indx = input('Index of new metabolite: '); % Choosing new metabolite to follow
        newMet = resultMets(indx,1);
        disp(' ')
        pos = find(resultMets(:,1) == newMet);
        if length(pos) == 1 % Metabolite used / created in one reaction
            rx = resultMets(pos,2);
        elseif length(pos) > 1 % Metabolite used / created in several reactions
            possibleRxns = resultMets(pos,2);
            for j = 1:length(possibleRxns)
                fprintf('%i:    %s\n\n', j, char(rxList(possibleRxns(j))))
            end
            indx = input('Choose index of correct reaction: ');
            rx = possibleRxns(indx);
        end
        if strcmp(dir, 'p') % Adding reaction to result (metabolite is product)
            if nargin == 1
                result(n,1) = cellstr(sprintf('%s ---> %s (rx %i: %s)', char(model.metNames(newMet)), char(model.metNames(metNr)), rx, char(model.rxnNames(rx))));
            elseif nargin == 2
                result(n,1) = cellstr(sprintf('%s -> (%8.3f) -> %s (rx %i: %s)', char(model.metNames(newMet)), FBAresult.x(rx), char(model.metNames(metNr)), rx, char(model.rxnNames(rx))));
            end
            n = n + 1;
        elseif strcmp(dir, 's') % Adding reaction to result (metabolite is substrate)
            if nargin == 1
                result(n,1) = cellstr(sprintf('%s ---> %s (rx %i: %s)', char(model.metNames(metNr)), char(model.metNames(newMet)), rx, char(model.rxnNames(rx))));
            elseif nargin == 2
                result(n,1) = cellstr(sprintf('%s -> (%8.3f) -> %s (rx %i: %s)', char(model.metNames(metNr)), FBAresult.x(rx), char(model.metNames(newMet)), rx, char(model.rxnNames(rx))));
            end
            n = n + 1;
        end
        metNr = newMet;
        lastRx = rx;
        clear rxns resultMets k otherMets pos newMet rx i j
        fprintf('\nReactions so far:\n\n')
        disp(char(result))
        disp(' ')
        go = input('Keep going? (yes = 1, no = 0): ');
        disp(' ')
    else
        disp('--- End of path ---')
        disp(' ')
        if strcmp(dir, 'p')
            result(n,1) = cellstr(sprintf('-| %s', char(model.metNames(metNr))));
        elseif strcmp(dir, 's')
            result(n,1) = cellstr(sprintf('%s -|', char(model.metNames(metNr))));
        end
        go = 0;
    end
    fprintf('-----\n\n')
end
disp(result)