newmapstack = reshape(mapStack,100,100,1000);
figure;
for i = 1:1000
    imagesc(newmapstack(:,:,i));axis image
    pause()
end


clearvars
close all
load('/Users/Garima/Documents/8980DM/project/dataset/Garima2/u1h99v99_1.mat')
n= 'u';
r= 1;
%mapStack= reshape(mapStack,100,100,1000);
mapStack= mapStack(:,1:r*100);

static_inds = sum(mapStack==1,2)==r*100 | sum(mapStack==2,2)==r*100;
R = mapStack(static_inds==0,:);
K = 16;
GT= GT(:,1:r*100);
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

pname=sprintf('final_data_results/%sP_%d_rep_%d.mat',n,K,r);
save(pname,'P');

qname=sprintf('final_data_results/%sQ_%d_rep_%d.mat',n,K,r);
save(qname,'Q');

%eR= P(:,7:16)*Q(7:16,:);
eR= P*Q;
T= findThresh(P,Q,remain_idx,R);
idx= cat(1,miss_idx,remain_idx);
T2= findThresh(P,Q,idx,R);
sameT= 0;
if T==T2
    sameT= 1;
end

%T= 1.5;
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
res= P*Q;
t = res(eRc==3 | eRc==4);
figure;hist(t);title('Histogram of Error Pixels after Thresholding');
xlabel('Predicted Values') ;
impute_res = zeros(size(data));
impute_res(data==1 & R==1) = 1;
impute_res(data==2 & R==2) = 2;
impute_res(data==1 & R==2) = 3;
impute_res(data==2 & R==1) = 4;

sprintf('one error: %d',size(eRc(eRc==3),1))
sprintf('two error: %d',size(eRc(eRc==4),1))
sprintf('impute one error: %d',size(impute_res(impute_res==3),1))
sprintf('impute two error: %d',size(impute_res(impute_res==4),1))

fmat2 = mapStack;
fmat2(static_inds==0,:) = eR;

figure;
for i = 1:r*100
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

figure,plotyy(1:1000,sum(mapStack==1,1),1:1000,sum(fmat2~=GT,1));



