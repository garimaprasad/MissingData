i
% load('/Users/Garima/Documents/8980DM/project/final_data_results/lP_16.mat')
% load('/Users/Garima/Documents/8980DM/project/final_data_results/lQ_16.mat')
% load('/Users/Garima/Documents/8980DM/project/dataset/Garima2/l1h99v99_1.mat')
% static_inds = sum(mapStack==1,2)==100 | sum(mapStack==2,2)==100;
% eR= P*Q;
% T= 1.5;
% eR(eR<=T)=1;
% eR(eR>T)=2;
% lMat = mapStack;
% lMat(static_inds==0,:) = eR;
% save('dataset/lMat.mat','lMat');


% bMap =load('/Users/Garima/Documents/8980DM/project/dataset/bMap.mat');
% bMap=bMap.maps;
% load('/Users/Garima/Documents/8980DM/project/dataset/bMat.mat');
% 
% rMap =load('/Users/Garima/Documents/8980DM/project/dataset/rMap.mat');
% rMap = rMap.maps;
% load('/Users/Garima/Documents/8980DM/project/dataset/rMat.mat');
% 
% uMap =load('/Users/Garima/Documents/8980DM/project/dataset/uMap.mat');
% uMap = uMap.mapu;
% load('/Users/Garima/Documents/8980DM/project/dataset/uMat.mat');
% 
% lMap =load('/Users/Garima/Documents/8980DM/project/dataset/lMap.mat');
% lMap = lMap.maps;
% load('/Users/Garima/Documents/8980DM/project/dataset/lMat.mat');
% 
% sMap =load('/Users/Garima/Documents/8980DM/project/dataset/sMap.mat');
% sMap = sMap.maps;
% load('/Users/Garima/Documents/8980DM/project/dataset/sMat.mat');


figure;
for i = 1:100
    subplot(3,2,1)
    imagesc(reshape(rMap(:,i),100,100));axis image
    title('Random Missing Data');
    subplot(3,2,2)
    t2 = rMat(:,i);
    t3 = GT(:,i);
    tf = t2;
    tf(t2==1 & t3==2) = 4;
    tf(t2==2 & t3==1) = 3;
    tf(1:4) = 1:4;
    imagesc(reshape(tf,100,100));axis image
    title('Predicted Data');
    
    subplot(3,2,3)
    imagesc(reshape(lMap(:,i),100,100));axis image
    title('Location Missing Data');
    
    subplot(3,2,4)
    t2 = lMat(:,i);
    t3 = GT(:,i);
    tf = t2;
    tf(t2==1 & t3==2) = 4;
    tf(t2==2 & t3==1) = 3;
    tf(1:4) = 1:4;
    imagesc(reshape(tf,100,100));axis image
    title('Predicted Data');
    
    subplot(3,2,5)
    imagesc(reshape(sMap(:,i),100,100));axis image
    title('Spatial Missing Data');
     
    subplot(3,2,6)
    t2 = sMat(:,i);
    t3 = GT(:,i);
    tf = t2;
    tf(t2==1 & t3==2) = 4;
    tf(t2==2 & t3==1) = 3;
    tf(1:4) = 1:4;
    imagesc(reshape(tf,100,100));axis image
    title('Predicted Data');
    
    pause()
end

figure;
for i= 1:100
    subplot(2,2,1)
    imagesc(reshape(uMap(:,i),100,100));axis image
    title('Spatio-Temporal Missing Data');
     
    subplot(2,2,2)
    t2 = uMat(:,i);
    t3 = GT(:,i);
    tf = t2;
    tf(t2==1 & t3==2) = 4;
    tf(t2==2 & t3==1) = 3;
    tf(1:4) = 1:4;
    imagesc(reshape(tf,100,100));axis image
    title('Predicted Data');
    
    subplot(2,2,3)
    imagesc(reshape(bMap(:,i),100,100));axis image
     title('Boundary Missing Data');
    
    subplot(2,2,4)
    t2 = bMat(:,i);
    t3 = GT(:,i);
    tf = t2;
    tf(t2==1 & t3==2) = 4;
    tf(t2==2 & t3==1) = 3;
    tf(1:4) = 1:4;
    imagesc(reshape(tf,100,100));axis image
    title('Predicted Data');
    pause()
end
