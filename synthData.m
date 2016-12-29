function [data] = synthData(rows,cols)
%tic;
data = ones(rows,cols);
%1 is water and 2 is land
for col = 1: cols
    k = round(rows*rand(1));
    data(1:k,col) = 2;
end
save('data.mat','data');
%time=toc;
% disp(time);
end