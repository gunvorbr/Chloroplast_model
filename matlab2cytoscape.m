function [] = matlab2cytoscape(model, filename, FBAresult)

disp(' ')

% Turning matlab model into csv-file which can be imported to cytoscape
%   by Gunvor Røkke, NTNU, 2021

time = clock; % Looking at the clock :-D
month = num2str(time(2));
if length(month) == 1
    month = sprintf('0%s', month);
end
filename = sprintf('%i_%s_%i_%s', time(1), month, time(3), filename); % Including time in file name
if isempty(strfind(filename, '.csv')) % Ensuring that filename is on correct format
    filename = sprintf('%s.csv', filename);
end

if nargin == 2
    network = cell(1,3);
elseif nargin == 3
    network = cell(1,4);
end
comp = input('Inducate letter of compartment, if you don''t want a compartment marked, press 0: ', 's');

disp(' ')
disp('Please be patient :-) This might take 20 seconds or so...')
disp(' ')

n = 1;

for i = 1:size(model.S,1) % Looping through metabolites
    for j = 1:size(model.S,2) % Looping through reactions
        if nonzeros(model.S(i,j)) < 0
            network(n,1) = model.metNames(i); % From metabolite (if met is used in rx)
            network(n,2) = model.rxns(j); % To rx
            network(n,3) = num2cell(model.rev(j)); % Reversibility
            if nargin == 3
                network(n,4) = num2cell(FBAresult.x(j));
            end
%            if isempty(strfind(comp, '0')) % Gathering compartment information about node interaction
%                metName = char(model.mets(i));
%                if strcmp(metName(end), ']')
%                    if strcmp(metName(end-1), comp)
%                        network(n,size(network,2)) = num2cell(1);
%                    else
%                        network(n,size(network,2)) = num2cell(0);
%                    end
%                elseif strcmp(metName(end-1), '_')
%                    if strcmp(metName(end), comp)
%                        network(n,size(network,2)) = num2cell(1);
%                    else
%                        network(n,size(network,2)) = num2cell(0);
%                    end
%                end
%            end
            n = n + 1;
        elseif nonzeros(model.S(i,j)) > 0
            network(n,1) = model.rxns(j); % From reaction (if met is produced)
            network(n,2) = model.metNames(i); % To met
            network(n,3) = num2cell(model.rev(j)); % Reversibility
            if nargin == 3
                network(n,4) = num2cell(FBAresult.x(j));
            end
%            if isempty(strfind(comp, '0')) % Gathering compartment information about node interaction
%                metName = char(model.mets(i));
%               if strcmp(metName(end), ']')
%                    if strcmp(metName(end-1), comp)
%                        network(n,size(network,2)) = num2cell(1);
%                    else
%                        network(n,size(network,2)) = num2cell(0);
%                    end
%                elseif strcmp(metName(end-1), '_')
%                    if strcmp(metName(end), comp)
%                        network(n,size(network,2)) = num2cell(1);
%                    else
%                        network(n,size(network,2)) = num2cell(0);
%                   end
%               end
%            end
            n = n + 1;
        end
    end
end
if strcmp(comp, '0') == 0 % Generating node attribute table
    compNodes = cell(length(model.mets),2);
    for i = 1:length(model.mets)
        compNodes(i,1) = model.metNames(i);
        metName = char(model.mets(i));
        if strcmp(metName(end), ']')
            if strcmp(metName(end-1), comp)
                compNodes(i,2) = num2cell(1);
            else
                compNodes(i,2) = num2cell(0);
            end
        elseif strcmp(metName(end-1), '_')
            if strcmp(metName(end), comp)
                compNodes(i,2) = num2cell(1);
            else
                compNodes(i,2) = num2cell(0);
            end
        end
    end
end

cd CytoscapeModels

cell2csv(filename, network, ';') % Saving network

if strcmp(comp, '0') == 0
    nodeFileName = sprintf('%s_compInfo%s', filename(1:end-4), filename(end-3:end));
    cell2csv(nodeFileName, compNodes, ';') % Saving node attribute
end

cd ..