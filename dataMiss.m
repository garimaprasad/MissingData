function[ data,data_miss]=dataMiss(replicate,noise)
%ground truth data
data= load('cf_gt_data.mat');
data= data.data;
%replicate data
x = data;
for i = 1:4
 x = cat(2,x,data);
end
data = x;
save('cf_gt_data_rep.mat','data');
rows=size(data,1);
cols=size(data,2);

%create missing data
%0 is missing data
%num_miss is of the order of the rows. This can be tuned to change the 
%level of sparseness and missing data.
data_miss= data;
for col = 1: cols
    num_miss = floor(rows*rand(1)*noise);
    temp = randperm(rows);   
    data_miss(temp(1:num_miss),col) = 0;
end

save('cf_Miss_data_replicate.mat','data_miss');
 %figure,imagesc(data_miss);title(sum(sum(data_miss==0))/(rows*cols))