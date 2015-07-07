function [datacell] = bootgrandaverage(STATS,condfiles_subs,numconds,nboot,numpnts,cond_bootvect,bootind,design,varargin)

%
% takes a resample of single subjects (their surrogate x time arrays), stacks them
% along 3rd dimension and takes the trimmed grand average. Different cases
% handle the different designs. Cases really just control how resampling
% indices change subtley with different designs.
%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if acceptable option was given
designopts={'bw','ww','bb', 'b', 'w'};
if ~any(strcmp(design, designopts));
    error('design input must be ''bw'', ''ww'', ''bb'', ''b'', or ''w''')
end

% check varargin length agains design length
if length(design)~=length(varargin);
    error(['If you specify a 2-way design (e.g., ''ww'') you must input jlvls and klvls as separate  trailing arguments' ...
        'For example, for a 2x4 design the last two inputs of the function would simply be, 2,4' ...
        'If you specify a 1-way design (e.g., ''w'') you must input only jlvls as the trailing agument' ...
        'For example, for a 1-way design with 3 levels the last input would be 3']);
end




if any(strcmp({'chanclust' 'gfa'},STATS.measure));
    
    %preallocating
    datacell=cell(1,size(condfiles_subs,2));
    
    for i=1:size(condfiles_subs,2);
        datacell{i}=zeros(nboot,numpnts);
    end
    
    switch design
        
        case 'bw'
            
            q=1;
            
            for j=1:jlvls;
                for k=1:klvls;
                    
                    [rowsubs ~]=size(condfiles_subs{1,q}(:,1));
                    subdata_gather=zeros(nboot,numpnts,rowsubs);
                    
                    for h=1:rowsubs;
                        tmpdata=load(condfiles_subs{1,q}{cond_bootvect.(['j',num2str(j)]).Data.dat(bootind,h)});
                        subdata_gather(:,:,h)=tmpdata.data;
                    end
                    
                    datacell{1,q}=trimmean(subdata_gather,40,3);
                    
                    q=q+1;
                end
                
            end
            
            clear subdata_gather tmpdata
            
        case {'ww','w'}
            
            for q=1:numconds;
                [rowsubs ~]=size(condfiles_subs{1,q}(:,1));
                subdata_gather=zeros(nboot,numpnts,rowsubs);
                
                for h=1:rowsubs;
                    tmpdata=load(condfiles_subs{1,q}{cond_bootvect.j1.Data.dat(bootind,h)});
                    subdata_gather(:,:,h)=tmpdata.data;
                end
                
                datacell{1,q}=trimmean(subdata_gather,40,3);
                
            end
            
            clear subdata_gather tmpdata
            
        case {'bb','b'}
            
            for j=1:numconds;
                
                [rowsubs ~]=size(condfiles_subs{1,j}(:,1));
                subdata_gather=zeros(nboot,numpnts,rowsubs);
                
                for h=1:rowsubs;
                    tmpdata=load(condfiles_subs{1,j}{cond_bootvect.(['j',num2str(j)]).Data.dat(bootind,h)});
                    subdata_gather(:,:,h)=tmpdata.data;
                end
                
                datacell{1,j}=trimmean(subdata_gather,40,3);
            end
            
            clear subdata_gather tmpdata
    end
    
elseif any(strcmp({'ersp' 'itc'},STATS.measure));
    
    % kill file
    warning off
    for i=1:length(STATS.condnames);
        delete(['mont_groupboots_',STATS.savestring, '_', STATS.condnames{i},'.map']);
    end
    warning on
    
    switch design
        
        case 'bw'
            
            q=1;
            for j=1:jlvls;
                for k=1:klvls;
                    [rowsubs ~]=size(condfiles_subs{1,q}(:,1));
                    subdata_gather=zeros(STATS.nboot,numpnts,rowsubs);
                    
                    for freqcurrent=1:STATS.freqbins;
                        
                        for h=1:rowsubs;
                            
                            
                            % memory map load
                            datamap=mapread(condfiles_subs{1,q}{cond_bootvect.(['j',num2str(j)]).Data.dat(bootind,h)},'dat','datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
                            subdata_gather(:,:,h)=squeeze(datamap.Data.dat(freqcurrent,:,:))';
                            
                        end
                        
                        mapwrite(trimmean(subdata_gather,40,3),['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'],'datsize',[STATS.nboot STATS.timesout STATS.freqbins]);
                        clear subdata_gather
                    end
                    
                    % reording array so taht everything maps correctly.
                    maptemp=mapread(['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'], 'dat');
                    mapwrite(permute(maptemp.Data.dat,[3 2 1]),['mont_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'],'datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
                    delete(['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map']);
                    
                    datacell{1,q}=mapread(['mont_groupboots_',STATS.savestring, '_', STATS.condnames{q},'.map'], 'dat');
                    q=q+1;
                end
            end
            
            
        case {'ww','w'}
            
            for q=1:numconds;
                [rowsubs ~]=size(condfiles_subs{1,q}(:,1));
                subdata_gather=zeros(STATS.nboot,numpnts,rowsubs);
                
                for freqcurrent=1:STATS.freqbins;
                    
                    for h=1:rowsubs;
                        
                        
                        % memory map load
                        datamap=mapread(condfiles_subs{1,q}{cond_bootvect.j1.Data.dat(bootind,h)},'dat');
                        subdata_gather(:,:,h)=squeeze(datamap.Data.dat(freqcurrent,:,:))';
                        
                    end
                    
                    mapwrite(trimmean(subdata_gather,40,3),['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'],'datsize',[STATS.nboot STATS.timesout STATS.freqbins]);
                    clear subdata_gather
                end
                
                % reording array so taht everything maps correctly.
                maptemp=mapread(['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'], 'dat');
                mapwrite(permute(maptemp.Data.dat,[3 2 1]),['mont_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'],'datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
                delete(['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map']);
                
                % filling data cell with objects
                datacell{1,q}=mapread(['mont_groupboots_',STATS.savestring, '_', STATS.condnames{q},'.map'], 'dat');
                
            end
            
            
            
        case {'bb','b'}
            
            for q=1:numconds;
                [rowsubs ~]=size(condfiles_subs{1,q}(:,1));
                subdata_gather=zeros(STATS.nboot,numpnts,rowsubs);
                
                for freqcurrent=1:STATS.freqbins;
                    
                    for h=1:rowsubs;
                        
                        
                        % memory map load
                        datamap=mapread(condfiles_subs{1,q}{cond_bootvect.(['j',num2str(q)]).Data.dat(bootind,h)},'dat','datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
                        subdata_gather(:,:,h)=squeeze(datamap.Data.dat(freqcurrent,:,:))';
                        
                    end
                    
                    mapwrite(trimmean(subdata_gather,40,3),['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'],'datsize',[STATS.nboot STATS.timesout STATS.freqbins]);
                    
                end
                
                maptemp=mapread(['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'], 'dat');
                mapwrite(permute(maptemp.Data.dat,[3 2 1]),['mont_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map'],'datsize',[STATS.freqbins STATS.timesout STATS.nboot]);
                delete(['temp_groupboots_',STATS.savestring,'_',STATS.condnames{q},'.map']);
                
                datacell{1,q}=mapread(['mont_groupboots_',STATS.savestring, '_', STATS.condnames{q},'.map'], 'dat');
                
            end
    end
    
    
    
end

end

