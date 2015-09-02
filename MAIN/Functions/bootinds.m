function [rowfile cond_bootvect tmpfname] = bootinds(condfiles_subs,nsamp,design,varargin)

%Creates bootstrap indices for given designs. Index arrays are memory mapped
%so as to not take up RAM.



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

switch design
    
    case 'bw'
        % create all random sampling indices
        q=1;
        k=1;
        for j=1:jlvls; %between lvls
            
            [rowfile(k) colfile]=size(condfiles_subs{1,q}(:,1));
            
            %%%%%% edit on may 8th/15
            [~,tmpfname_tmp]=fileparts(tempname);
            tmpfname{k}=tmpfname_tmp;
            %%%%%%
            
            %%%%%% edit on may 8th/15
            %fidm=mapwrite(randi(rowfile(k),nsamp,rowfile(k)),['cond_bootvect',num2str(k),'.map'],'datsize',[nsamp rowfile(k)]);
            fidm=mapwrite(randi(rowfile(k),nsamp,rowfile(k)),[tmpfname{k},'.map'],'datsize',[nsamp rowfile(k)]);
            
            q=q+klvls;
            k=k+1;
        end
        
        for i=1:jlvls
            cond_bootvect.(['j',num2str(i)])=mapread([tmpfname{k},'.map'],'dat','datsize',[nsamp rowfile(i)]);
        end
        
        
    case {'ww','w'}
        
        % this assumes that for within subjects designs, the number of subjects is equal in all cells.
        [rowfile colfile]=size(condfiles_subs{1,1}(:,1));
        
        [~,tmpfname]=fileparts(tempname); % unique file name
        %fidm=mapwrite(randi(rowfile,nsamp,rowfile),'cond_bootvect1.map','datsize',[nsamp rowfile]);
        fidm=mapwrite(randi(rowfile,nsamp,rowfile),[tmpfname, '.map'],'datsize',[nsamp rowfile]);
        
        %cond_bootvect.j1=mapread('cond_bootvect1.map','dat','datsize',[nsamp rowfile]);
        cond_bootvect.j1=mapread([tmpfname, '.map'],'dat','datsize',[nsamp rowfile]);
        
    case {'bb','b'}
        
        % create different random sampling indices for each condition
        for j=1:length(condfiles_subs);
            %for j=1:numconds;
            
            [rowfile(j) colfile]=size(condfiles_subs{1,j}(:,1));
            
             %%%%%% edit on may 8th/15
            [~,tmpfname_tmp]=fileparts(tempname);
            tmpfname{j}=tmpfname_tmp;
            %%%%%%
            
            fidm=mapwrite(randi(rowfile(j),nsamp,rowfile(j)),[tmpfname{j},'.map'],'datsize',[nsamp rowfile(j)]);
            
        end
        
        for i=1:length(condfiles_subs);
            %i=1:numconds
            cond_bootvect.(['j',num2str(i)])=mapread([tmpfname{i},'.map'],'dat','datsize',[nsamp rowfile(i)]);
        end
        
    
end
end
