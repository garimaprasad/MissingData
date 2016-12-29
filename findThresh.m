function T= findThresh(P,Q,remain_idx,data)
% nbins=200;
% figure;
% histogram(eR-1,nbins);
% [counts,x] = imhist(eR-1,nbins);
% T = graythresh(counts);
min_err=length(remain_idx);
T=0;
for t = 1:0.1:2
eR = P*Q;
eR(eR<=t)=1;
eR(eR>t)=2;
res = zeros(size(data));
res(data(remain_idx)==1 & eR(remain_idx)==1) = 1;
res(data(remain_idx)==2 & eR(remain_idx)==2) = 2;
res(data(remain_idx)==1 & eR(remain_idx)==2) = 3;
res(data(remain_idx)==2 & eR(remain_idx)==1) = 4;
error = size(res(res==4),1)+size(res(res==3),1);
if error <= min_err
    min_err = error;
    T = t;
end
end
sprintf('Threshold: %d',T);
end