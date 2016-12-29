function [rem_1,rem_2,miss_one,miss_two,n_per,imp_1_err,imp_2_err,rem_1_err,rem_2_err,error_impute]=test_cf(data,noise,K)
rows=size(data,1);
cols=size(data,2);

% %create missing data
% %0 is missing data
% %num_miss is of the order of the rows. This can be tuned to change the 
% %level of sparseness and missing data.
data_miss= data;
for col = 1: cols
    num_miss = floor(rows*rand(1)*noise);
    temp = randperm(rows);   
    data_miss(temp(1:num_miss),col) = 0;
end
n_per=sum(sum(data_miss==0))*100/(rows*cols);
%figure,imagesc(data_miss);title('noise percent: %d',sum(sum(data_miss==0))*100/(rows*cols))

R = data_miss;
miss_idx = find(R==0);
miss_one= sum(data(miss_idx)==1);
miss_two= sum(data(miss_idx)==2);
rem_1=sum(sum(data_miss==1));
rem_2=sum(sum(data_miss==2));

R= transpose(R);
%initialize P and Q
[Pn,Qn] = nnmf(R,K);
Qn = Qn';
[P,Q] = matrix_factorization(R, Pn, Qn, K);
eR= P*Q;
eR=transpose(eR);

%figure; hist(eR,100)
% T= findThresh(eR);
% disp(T);
T=1.5;
eR(eR<=T)=1;
eR(eR>T)=2;

data_miss(miss_idx)= eR(miss_idx);

eRc = zeros(size(data));
eRc(data==1 & eR==1) = 1;
eRc(data==2 & eR==2) = 2;
eRc(data==1 & eR==2) = 3;
eRc(data==2 & eR==1) = 4;

impute_res = zeros(size(data));
impute_res(data==1 & data_miss==1) = 1;
impute_res(data==2 & data_miss==2) = 2;
impute_res(data==1 & data_miss==2) = 3;
impute_res(data==2 & data_miss==1) = 4;

imp_1_err=size(impute_res(impute_res==3),1);
imp_2_err=size(impute_res(impute_res==4),1);
rem_1_err=size(eRc(eRc==3),1)-imp_1_err;
rem_2_err=size(eRc(eRc==4),1)-imp_2_err;
error_impute=size(find((data(miss_idx)-eR(miss_idx)) ~=0),1);
end

function [P,Q] = matrix_factorization(R, P, Q, K)
steps=1000;
alpha=0.0002;
beta=0.02;
Q= transpose(Q);
eprev=0;
for step = 1: steps
  
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
     e = sum(sum(R-P*Q));
     e=e+(beta/2) * (power(norm(P),2) + power(norm(transpose(Q)),2));
     
     if abs(e-eprev) < 0.0001
        break;
     end
     eprev=e;
end
end