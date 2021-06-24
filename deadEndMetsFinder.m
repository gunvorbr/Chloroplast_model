function deadEndList = deadEndMetsFinder(model)

% Identifies metabolites in a model that are only produced or only consumed
%   by Gunvor Røkke, NTNU, 2021

deadEndList = double.empty(0,1);
n = 1;

for i = 1:length(model.mets)
    rxns = find(model.S(i,:) ~= 0);
    if length(rxns) == 1
        deadEndList(n,1) = i;
        n = n + 1;
    elseif length(rxns) > 1
        if isequal(model.rev(rxns), zeros(length(rxns),1)) % All rxns are irreversible
            Scoeffs = nonzeros(model.S(i,rxns));
            if (length(find(Scoeffs > 0)) == length(rxns)) || (length(find(Scoeffs < 0)) == length(rxns))
                deadEndList(n,1) = i;
                n = n + 1;
            end
        end
    end
end