function [windata] = winsorize(data,tr)
% winsorize columns of input matrix or vector (data)
% 20% winsorising (.2) is usually recommended, however, you must specify this in the input

% winsorize data
[rows cols]=size(data);
winind=round(tr*rows);
windata=data;
data_sort=sort(data);

for i=1:cols;
    winlow=data(:,i)<=data_sort(winind+1,i);
    winhi=data(:,i)>=data_sort(end-winind,i);
    windata(winlow,i)=data_sort(winind+1,i);
    windata(winhi,i)=data_sort(end-winind,i);
    clear winlow winhi
end

end

