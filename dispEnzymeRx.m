function [] = dispEnzymeRx(ECnum)

% Searches through the KEGG database, and displays the reactions of one particular enzyme
%   by Gunvor Røkke, NTNU, 2021

res = urlread(strcat('http://rest.kegg.jp/get/', ECnum));

% Checking if new EC-number has been assigned

if isempty(strfind(res, 'ALL_REAC')) && ~isempty(strfind(res, 'Transferred to'))
    newEC = strfind(res, 'Transferred to');
    if length(newEC) == 1
        newEC = newEC + 15;
    else
        newEC = newEC(1) + 15;
    end
    pointMarker = zeros(3,1);
    pointFinder = newEC;
    n = 1;
    while ~isequal(pointMarker, ones(3,1))
        pointFinder = pointFinder + 1;
        if strcmp(res(pointFinder), '.')
            pointMarker(n,1) = 1;
            n = n + 1;
        end
    end
    clear pointMarker n
    newECend = pointFinder + 1;
    clear pointFinder
    while ~isnan(str2double(res(newECend)))
        newECend = newECend + 1;
    end
    newECend = newECend - 1;
    res = urlread(strcat('http://rest.kegg.jp/get/', res(newEC:newECend)));
end

% Searching for reactions 

rxStart = strfind(res, 'ALL_REAC');
if isempty(rxStart)
    fprintf('\nNo reactions found\n\n')
else
    rxStart = rxStart + 12;
    rxStopp = rxStart;
    go = 1;
    while go == 1
        if strcmp(res(rxStopp), sprintf('\n')) && ~strcmp(res(rxStopp+1), ' ')
            go = 0;
        else
            rxStopp = rxStopp + 1;
        end
    end
    rxStopp = rxStopp - 1;
    rawRxns = res(rxStart:rxStopp);
    rawRxns = strrep(rawRxns, ';', '');
    rawRxns = strsplit(rawRxns);
    if ~isempty(smatch(rawRxns, '(other)'))
        entry = smatch(rawRxns, '(other)');
        rawRxns(entry) = [];
    end
    rxIndicator = smatch(rawRxns, 'R');
    rxns = rawRxns(rxIndicator);
    for i = 1:length(rxns)
        fprintf('\nRx %i:\n', i)
        dispKEGGrx(char(rxns(i)));
        disp(' ')
    end
end