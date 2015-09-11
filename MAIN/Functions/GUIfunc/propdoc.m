function [proptext] = propdoc(func,startstr,endstr)

a=help(func);
b=strfind(a,startstr);
c=a(b+length(startstr):end);
d=strfind(c,endstr);
proptext=c(1:d(1)-1);

end

