clearvars
close all
load('/Users/Garima/Documents/8980DM/project/dataset/Garima2/u2h99v99_1.mat')
noise= 'u2';
static_inds = sum(mapStack==1,2)==100 | sum(mapStack==2,2)==100;
R = mapStack(static_inds==0,:);
load('/Users/Garima/Documents/8980DM/project/dataset/BaselineMethods/u2h99v99_1.mat','fmapStack_withoutSnow');      
K = 16;
data= GT(static_inds==0,:);
miss_idx = find(R==0);
remain_idx= find(R>0);
eR = fmapStack_withoutSnow(static_inds==0,:);
eR(eR==4)= 2;
R(miss_idx)= eR(miss_idx);

eRc = zeros(size(data));
eRc(data==1 & eR==1) = 1;
eRc(data==2 & eR==2) = 2;
eRc(data==1 & eR==2) = 3;
eRc(data==2 & eR==1) = 4;


impute_res = zeros(size(data));
impute_res(data(miss_idx)==1 & eR(miss_idx)==1) = 1;
impute_res(data(miss_idx)==2 & eR(miss_idx)==2) = 2;
impute_res(data(miss_idx)==1 & eR(miss_idx)==2) = 3;
impute_res(data(miss_idx)==2 & eR(miss_idx)==1) = 4;

sprintf('one error: %d',size(eRc(eRc==3),1))
sprintf('two error: %d',size(eRc(eRc==4),1))
sprintf('impute one error: %d',size(impute_res(impute_res==3),1))
sprintf('impute two error: %d',size(impute_res(impute_res==4),1))

fmat2 = mapStack;
fmat2(static_inds==0,:) = eR;

figure;
for i = 1:100
    subplot(1,3,1)
    imagesc(reshape(mapStack(:,i),100,100));axis image
    
    subplot(1,3,2)
    t1 = mapStack(:,i);
    t2 = fmat2(:,i);
    t3 = GT(:,i);
    tf = t2;
    tf(t2==1 & t3==2) = 4;
    tf(t2==2 & t3==1) = 3;
    tf(1:4) = 1:4;
    imagesc(reshape(tf,100,100));axis image
    
    subplot(1,3,3)
    imagesc(reshape(GT(:,i),100,100));axis image
    pause()
end

figure,plotyy(1:100,sum(mapStack==1,1),1:100,sum(fmat2==1&GT==2,1));
title('Time Series of Water Label');
figure,plotyy(1:100,sum(mapStack==2,1),1:100,sum(fmat2==2&GT==1,1));
title('Time Series of Water Label');


[~,ix] = sort(sum(data==1,2),'ascend');
figure,imagesc(eR(ix,:))
newR = mapStack(static_inds==0,:);
figure,imagesc(newR(ix,:))


% clearvars
% load('/Users/Garima/Documents/8980DM/project/dataset/BaselineMethods/u2h99v99_1.mat');
% figure;imagesc(reshape(dyn_inds,100,100))

