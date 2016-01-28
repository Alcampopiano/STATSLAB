function [hist_str]=statslab_history(varargin)

% determine calling function
stDebug = dbstack;
funcname=stDebug(2).name;

opts_pre=varargin(1:end-1); % get everything but the varargin at the end
opts_var=varargin(end); % get varargin

hist_str=['[STATS]=', funcname, '('];
outstr_pre=vararg2str_statslab(opts_pre);

% separate vararginou
outstr_var=vararg2str_statslab(opts_var{1});

% build string
hist_str=[hist_str outstr_pre, ',', outstr_var, ');'];

