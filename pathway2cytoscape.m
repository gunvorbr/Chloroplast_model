function [] = pathway2cytoscape(model, vector)

% Turning matlab model into csv-file for import in cytoscape
%   by Gunvor Røkke, NTNU, 2021

disp(' ')

time = clock; % Looking at the clock :-D
month = num2str(time(2));
if length(month) == 1
    month = sprintf('0%s', month);
end
day = num2str(time(3));
if length(day) == 1
    day = sprintf('0%s', day);
end
filename = input('Filename? ', 's');
filename = sprintf('%i_%s_%s_%s', time(1), month, day, filename); % Including time in file name
if isempty(strfind(filename, '.csv')) % Ensuring that filename is on correct format
    filename = sprintf('%s.csv', filename);
end
disp(' ')

% Checking format of vector and converting if necessary

if iscell(vector)
    newVector = zeros(length(vector),1);
    for i = 1:length(vector)
        newVector(i,1) = smatch(char(vector(i)), model.rxns, 'exact');
    end
    vector = newVector;
    vector = sort(vector);
    clear newVector
end

% Creating variables

network = cell(1,3);

n = 1;

for j = 1:length(vector) % Looping through reactions in vector
    for i = 1:size(model.S,1) % Looping through metabolites
        if nonzeros(model.S(i,vector(j))) < 0
            network(n,1) = model.metNames(i); % From metabolite (if met is used in rx)
            network(n,2) = model.rxns(vector(j)); % To rx
            network(n,3) = num2cell(model.rev(vector(j))); % Reversibility
            n = n + 1;
        elseif nonzeros(model.S(i,vector(j))) > 0
            network(n,1) = model.rxns(vector(j)); % From reaction (if met is produced)
            network(n,2) = model.metNames(i); % To met
            network(n,3) = num2cell(model.rev(vector(j))); % Reversibility
            n = n + 1;
        end
    end
end

cd CytoscapeModels

cell2csv(filename, network, ';') % Saving network

cd ..