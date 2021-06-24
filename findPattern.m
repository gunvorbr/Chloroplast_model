function result = findPattern(cellArray, pattern)

% Searches all rows of cell array to find out if a certain pattern is
% present.
%   by Gunvor Røkke, NTNU, 2021

result = cell(1,2);
n = 1;

for i = 1:length(cellArray)
    if strfind(lower(char(cellArray(i))), lower(char(pattern)))
        result(n,1) = num2cell(i);
        result(n,2) = cellArray(i);
        n = n + 1;
    end
end

