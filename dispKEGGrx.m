function [] = dispKEGGrx(KEGGID)

% Displays info about KEGG metabolites (search with KEGG number)
% Lists info as reactions
% If input is a reaction, displays reaction info
%   by Gunvor Røkke, NTNU, 2021

res = urlread(strcat('http://rest.kegg.jp/get/', KEGGID));

if strcmp(KEGGID(1), 'C')
    if ~isempty(res)
        start = strfind(res, 'REACTION') + 12;
        if strfind(res, 'PATHWAY')
            slutt = strfind(res, 'PATHWAY') - 2;
        elseif strfind(res, 'MODULE')
            slutt = strfind(res, 'MODULE') - 2;
        elseif strfind(res, 'ENZYME')
            slutt = strfind(res, 'ENZYME') - 2;
        elseif strfind(res, 'BRITE')
            slutt = strfind(res, 'BRITE') - 2;
        elseif strfind(res, 'DBLINKS')
            slutt = strfind(res, 'DBLINKS') - 2;
        end
        rxns = res(start:slutt);
        rxList = transpose(strsplit(rxns));
        for i = 1:length(rxList)
            res = urlread(strcat('http://rest.kegg.jp/get/', char(rxList(i))));
            start1 = strfind(res, 'DEFINITION') + 12;
            slutt1 = strfind(res, 'EQUATION') - 2;
            start2 = strfind(res, 'EQUATION') + 12;
            if strfind(res, 'COMMENT')
                slutt2 = strfind(res, 'COMMENT') - 2;
            elseif strfind(res, 'RPAIR')
                slutt2 = strfind(res, 'RPAIR') - 2;
            elseif strfind(res, 'ENZYME')
                slutt2 = strfind(res, 'ENZYME') - 2;
            elseif strfind(res, 'PATHWAY')
                slutt2 = strfind(res, 'PATHWAY') - 2;
            elseif strfind(res, 'MODULE')
                slutt2 = strfind(res, 'MODULE') - 2;
            elseif strfind(res, 'ORTHOLOGY')
                slutt2 = strfind(res, 'ORTHOLOGY') - 2;
            elseif strfind(res, 'RCLASS')
                slutt2 = strfind(res, 'RCLASS') - 2;
            elseif strfind(res, 'DBLINKS')
                slutt2 = strfind(res, 'DBLINKS') - 2;
            elseif strfind(res, '///')
                slutt2 = strfind(res, '///') - 2;
            end
            if strfind(res, 'ENZYME')
                ECstart = strfind(res, 'ENZYME') + 12;
                ECslutt = length(res);
                EClist = strsplit(res(ECstart:ECslutt));
                EC = 1;
            else
                EC = 0;
            end
            disp(' ')
            fprintf('Reaction %i, %s:', i, char(rxList(i)))
            disp(' ')
            disp(res(start1:slutt1))
            disp(res(start2:slutt2))
            if EC == 1
                disp(EClist(1))
            end
            disp(' ')
        end
    elseif isempty(res)
        warning('Uh oh...')
    end
elseif strcmp(KEGGID(1), 'R')
    if strfind(res, 'ENTRY')
        iStart = strfind(res, 'ENTRY') + 12;
        iStopp = iStart + 5;
        fprintf('\nReaction ID:         %s\n', res(iStart:iStopp))
    end
    if strfind(res, 'NAME')
        nStart = strfind(res, 'NAME') + 12;
        nStopp = nStart;
        while ~strcmp(res(nStopp), sprintf('\n'))
            nStopp = nStopp + 1;
        end
        nStopp = nStopp - 1;
        fprintf('Reaction name:       %s\n', res(nStart:nStopp))
    end
    if strfind(res, 'DEFINITION')
        eStart = strfind(res, 'DEFINITION') + 12;
        eStopp = eStart;
        while ~strcmp(res(eStopp), sprintf('\n'))
            eStopp = eStopp + 1;
        end
        eStopp = eStopp - 1;
        fprintf('Reaction equation:   %s\n', res(eStart:eStopp))
    end
    if strfind(res, 'EQUATION')
        kStart = strfind(res, 'EQUATION') + 12;
        kStopp = kStart;
        while ~strcmp(res(kStopp), sprintf('\n'))
            kStopp = kStopp + 1;
        end
        kStopp = kStopp - 1;
        fprintf('KEGG ID equation:    %s\n\n', res(kStart:kStopp))
    end
end    