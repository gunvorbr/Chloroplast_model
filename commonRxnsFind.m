function [] = commonRxnsFind(model)

% Script lets user imput two metabolite names, and searches
% through S matrix to see if there are reactions in the model that both
% metabolites participate in
%   by Gunvor Røkke, NTNU, 2021

rxList = ReactionNames(model);

n = 1;

met1 = input('Metabolite 1? ', 's');
met1List = findPattern(model.metNames, met1);
disp([transpose(num2cell(1:size(met1List,1))), met1List])
input1 = input('Choose number of correct metabolite: ');
disp(' ')
met1i = cell2mat(met1List(input1,1));


met2 = input('Metabolite 2? ', 's');
met2List = findPattern(model.metNames, met2);
disp([transpose(num2cell(1:size(met2List,1))), met2List])
input2 = input('Choose number of correct metabolite: ');
disp(' ')
met2i = cell2mat(met2List(input2,1));

if ~isempty(met1i) && ~isempty(met2i)
    for i = 1:size(model.S,2)
        if (model.S(met1i,i) ~= 0) && (model.S(met2i,i) ~= 0)
            rxNum(n,1) = i;
            n = n + 1;
        end
    end
    if exist('rxNum')
        for i = 1:length(rxNum)
            excRx(i,1) = rxList(rxNum(i));
        end
        disp([num2cell(rxNum), excRx])
    else
        disp('The two metabolites have no common reactions')
    end
end