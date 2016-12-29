function error_impute=cf(replicate,noise,K)
%ground truth data
% data= load('cf_gt_data_rep.mat');
% data= data.data;
% data_miss= load('cf_Miss_data_replicate.mat');
% data_miss= data_miss.data_miss;
%data= synthData(100,20);
data= load('data.mat');
 data= data.data;
%replicate data
% x = data;
% for i = 1:replicate
%  x = cat(2,x,data);
% end
% data = x;

rows=size(data,1);
cols=size(data,2);
% 
% %create missing data
% %0 is missing data
% %num_miss is of the order of the rows. This can be tuned to change the 
% %level of sparseness and missing data.
data_miss= data;
for col = 1: cols
    num_miss = floor(rows*rand(1)*0.4);
    temp = randperm(rows);   
    data_miss(temp(1:num_miss),col) = 0;
end
noiseper=sum(sum(data_miss==0))*100/(rows*cols)
figure,imagesc(data_miss);title(['Dat with percent: ',num2str(noiseper),' %'])

R = data_miss;
%R = data;
miss_idx = find(R==0);
remain_idx= find(R>0);
miss_one= sum(data(miss_idx)==1)
miss_two= sum(data(miss_idx)==2)
R= transpose(R);
N=size(R,1);
M= size(R,2);
P = rand(N,K);
Q = rand(M,K);
[Pn,Qn] = nnmf(R,K);
Qn = Qn';
[P,Q] = matrix_factorization(R, Pn, Qn, K);
eR= P*Q;
eR=transpose(eR);
%save('cf_eR_B4T_rep_slow.mat','eR');
%figure; hist(eR,100)
T= findThresh(P,Q,remain_idx,data);
% disp(T);
%T=1.5;
eR(eR<=T)=1;
eR(eR>T)=2;
%save('cf_eR_afterT_rep_slow.mat','eR');
data_miss(miss_idx)= eR(miss_idx);
%save('cf_result_rep_slow.mat','data_miss');
eRc = zeros(size(data));
eRc(data==1 & eR==1) = 1;
eRc(data==2 & eR==2) = 2;
eRc(data==1 & eR==2) = 3;
eRc(data==2 & eR==1) = 4;
t=[];
%histogram of error pixels
res= P*Q;
res=transpose(res);
t = res(eRc==3 | eRc==4);
figure;hist(t);


impute_res = zeros(size(data));
impute_res(data==1 & data_miss==1) = 1;
impute_res(data==2 & data_miss==2) = 2;
impute_res(data==1 & data_miss==2) = 3;
impute_res(data==2 & data_miss==1) = 4;

%save('cf_res_encoded_rep_slow.mat','eRc');
sprintf('one error: %d',size(eRc(eRc==3),1))
sprintf('two error: %d',size(eRc(eRc==4),1))
% sprintf('total error: %d',size(eRc(eRc==4),1)+size(eRc(eRc==3),1))
error_impute=size(find((data(miss_idx)-eR(miss_idx)) ~=0),1);
% sprintf('total missing values: %d',size(miss_idx,1))
%sprintf('total imputes is: %d ', size(miss_idx,1))
sprintf('impute one error: %d',size(impute_res(impute_res==3),1))
sprintf('impute two error: %d',size(impute_res(impute_res==4),1))
% sprintf('impute error is: %d ', error_impute)
 figure,imagesc(impute_res);title('Error Encoded impute result');
 figure,imagesc(data);title('original data');
 figure,imagesc(data_miss);title('Result');
 figure,imagesc(eRc);title('Error Encoded result');

end

function [P,Q] = matrix_factorization(R, P, Q, K)
Errors = [];

steps=10000;
alpha=0.0002;
beta=0.02;
Q= transpose(Q);
eprev=0;
for step = 1: steps
%      E = R - P*Q;
%      P = P + alpha * (2 * E * transpose(Q) - beta * P);
%      Q = Q + alpha * (2 * transpose(P)* E - beta * Q);    
    for i = 1: size(R,1)
        for j =1 : size(R,2)
            if R(i,j) > 0
                eij = R(i,j) - (P(i,:)*Q(:,j));
                for k =1:K
                    P(i,k) = P(i,k) + alpha * (2 * eij * Q(k,j) - beta * P(i,k));
                    Q(k,j) = Q(k,j) + alpha * (2 * eij * P(i,k) - beta * Q(k,j));
                end
            end
        end
    end
    e=0;
    for i = 1: size(R,1)
        for j =1 : size(R,2)
           if R(i,j) > 0
                e = e+ power((R(i,j) - (P(i,:)*Q(:,j))),2); 
                for k =1:K
                    e= e+(beta/2)* (power(P(i,k),2) + power(Q(k,j),2));
                end
           end
        end
    end
%      e = sum(sum(R-P*Q));
%      e=e+(beta/2) * (power(norm(P),2) + power(norm(transpose(Q)),2));
     Errors = [Errors e];
     if abs(e-eprev) < 0.0001
        break;
     end
     eprev=e;
end
step
%save('P_replicate_slow.mat','P');
%save('Q_replicate_slow.mat','Q');
figure,plot(Errors);
end