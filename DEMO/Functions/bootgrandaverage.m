function [datacell] = bootgrandaverage(condfiles_subs,numconds,nboot,numpnts,cond_bootvect,bootind,design,varargin)

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
end

