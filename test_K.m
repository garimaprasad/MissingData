clearvars
close all
load('/Users/Garima/Documents/8980DM/project/dataset/Garima2/u2h99v99_1.mat')
noise= 'u2';
static_inds = sum(mapStack==1,2)==100 | sum(mapStack==2,2)==100;
R = mapStack(static_inds==0,:);
% rn_inds = rand(size(R)>0.8 & R~=0);
% rn_inds= randi(numel(remain_idx), 0.2*numel(M));
% R(rn_inds) = 3 - R(rn_inds);
% newMapStack= mapstack;
% newMapStack(static_inds==0,:) = R;
K = 12;
data= GT(static_inds==0,:);
miss_idx = find(R==0);
remain_idx= find(R>0);
[Pn,Qn] = nnmf(R,K);
%fmat = Pn*Qn;
%fmat2 = mapStack;
%fmat2(static_inds==0,:) = fmat;
Qn = Qn';
steps=1000;

[P,Q] = matrix_factorization(R, Pn, Qn, K, steps);

pname=sprintf('final_data_results/%sP_%d.mat',noise,K);
save(pname,'P');

qname=sprintf('final_data_results/%sQ_%d.mat',noise,K);
save(qname,'Q');

eR= P*Q;

T= findThresh(P,Q,remain_idx,R);
% idx= cat(1,miss_idx,remain_idx);
% T2= findThresh(P,Q,idx,R);
% sameT= 0;
% if T==T2
%     sameT= 1;
% end
eR(eR<=T)=1;
eR(eR>T)=2;
R(miss_idx)= eR(miss_idx);

eRc = zeros(size(data));
eRc(data==1 & eR==1) = 1;
eRc(data==2 & eR==2) = 2;
eRc(data==1 & eR==2) = 3;
eRc(data==2 & eR==1) = 4;
t=[];
%histogram of error pixels

% res= P*Q;
% t = res(eRc==3 | eRc==4);
% mFigure = figure();
% mTextBox = uicontrol('style','text');
% mTextBox1 = uicontrol('style','text');
% position1= [1000 300 100 20];
% position= [1000 350 40 20];
% set(mTextBox,'String','K=16');
% set(mTextBox, 'Position', position);
% %set(mTextBox1, 'String', 'Spatio-temporal Holes', 'Position', position1);
% hist(t);

impute_res = zeros(size(data));
impute_res(data(miss_idx)==1 & eR(miss_idx)==1) = 1;
impute_res(data(miss_idx)==2 & eR(miss_idx)==2) = 2;
impute_res(data(miss_idx)==1 & eR(miss_idx)==2) = 3;
impute_res(data(miss_idx)==2 & eR(miss_idx)==1) = 4;

sprintf('one error: %d',size(eRc(eRc==3),1))
sprintf('two error: %d',size(eRc(eRc==4),1))
sprintf('impute one error: %d',size(impute_res(impute_res==3),1))
sprintf('impute two error: %d',size(impute_res(impute_res==4),1))
total_err= (size(eRc(eRc==3),1) +size(eRc(eRc==4),1))*100/numel(data)
imp_err= (size(impute_res(impute_res==3),1)+size(impute_res(impute_res==4),1))*100/numel(miss_idx)

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


