% Generates a similar list as the function printRxnFormula in COBRA
% toolbox, but uses metabolite names instead of IDs.

function [RxnEquNames] = ReactionNames(model)

RxnEquNames = cell(length(model.rxns),1);

for i = 1:length(model.rxns)
    if ~isempty(find(model.S(:,i) < 0)) % Reaction contains substrate(s)
        usedList = find(model.S(:,i) < 0);
        usedCoeff = zeros(length(usedList),1);
        for j = 1:length(usedCoeff)
            usedCoeff(j) = model.S(usedList(j),i);
        end
        clear j
        usedNames = cell(length(usedList),1);
        for j = 1:length(usedList)
            usedNames(j) = model.metNames(usedList(j));
        end
        clear j
    end
    if ~isempty(find(model.S(:,i) > 0)) % Reaction contains product(s)
        producedList = find(model.S(:,i) > 0);
        producedCoeff = zeros(length(producedList),1);
        for j = 1:length(producedList)
            producedCoeff(j) = model.S(producedList(j),i);
        end
        clear j
        producedNames = cell(length(producedList),1);
        for j = 1:length(producedList)
            producedNames(j) = model.metNames(producedList(j));
        end
        clear j
    end
    if exist('usedNames')
        usedText = sprintf('%g %s', (-1*usedCoeff(1)), char(usedNames(1)));
        for j = 2:length(usedList)
            usedText = sprintf('%s + %g %s', usedText, (-1*usedCoeff(j)), char(usedNames(j)));
        end
        clear j
    else
        usedText = ' ';
    end
    if exist('producedNames')
        producedText = sprintf('%g %s', producedCoeff(1), char(producedNames(1)));
        for j = 2:length(producedList)
            producedText = sprintf('%s + %g %s', producedText, producedCoeff(j), char(producedNames(j)));
        end
        clear j
    else
        producedText = ' ';
    end
    if (model.ub(i) > 0) && (model.lb(i) == 0)
        finalText = sprintf('%s -> %s', usedText, producedText);
    elseif (model.ub(i) == 0) && (model.lb(i) < 0)
        finalText = sprintf('%s <- %s', usedText, producedText);
    elseif (model.ub(i) > 0) && (model.lb(i) < 0)
        finalText = sprintf('%s <=> %s', usedText, producedText);
    elseif (model.ub(i) > 0) && (model.lb(i) > 0)
        finalText = sprintf('%s -> %s', usedText, producedText);
    elseif (model.ub(i) < 0) && (model.lb(i) < 0)
        finalText = sprintf('%s <- %s', usedText, producedText);
    elseif (model.ub(i) == 0) && (model.lb(i) == 0)
        finalText = sprintf('%s --- %s', usedText, producedText);
    end
    RxnEquNames(i) = cellstr(finalText);
    clear usedList usedCoeff usedNames producedList producedCoeff producedNames usedText producedText
end