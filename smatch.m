function hits = smatch(array, sTerm, exact)

% Works exactly as strmatch, but saves you of the annoying orange warning
% lines when you are using it in a script

%   by Gunvor Røkke, NTNU, 2021

hits = double.empty(0,1);
n = 1;

if nargin == 2 % Searching for pattern
    for i = 1:length(array)
        if strfind(lower(char(array(i))), lower(sTerm))
            hits(n,1) = i;
            n = n + 1;
        end
    end
elseif (nargin == 3) && strcmp(exact, 'exact') % Searching for exact match
    for i = 1:length(array)
        if strcmp(char(array(i)), sTerm)
            hits(n,1) = i;
            n = n + 1;
        end
    end
end