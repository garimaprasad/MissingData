function eRc = bpr_stochastic()

data= load('cf_gt_data_rep.mat');
data= data.data;
dmiss= load('cf_Miss_data_replicate.mat');
dmiss= dmiss.data_miss;

%make matrix 
users= size(dmiss,2);
items= size(dmiss,1);
figure,imagesc(dmiss);title(sum(sum(dmiss==0))/(items*users))
miss_idx = find(dmiss==0);


K= 3;
P = load('P_replicate_slow.mat');
P=P.P;
Q = load('Q_replicate_slow.mat');
Q=transpose(Q.Q);
last_loss=0;
lRate= 0.001;
reg=0.001;
x = P* transpose(Q);           

steps=0;                
while steps <= 5000
    loss=0;
    u= randi([1 users]);
    i=randi([1 items]);
    j=randi([1 items]);
        if i ~=j
            steps= steps+1;
            xuij = x(u,i)-x(u,j);
            loss = loss + log(1/(1+exp(-1*xuij)));
            sigmoid= 1/(1+exp(xuij));

            for k= 1:K
                puk= P(u,k);
                qik= Q(i,k);
                qjk= Q(j,k);
                P(u,k) = P(u,k) - lRate * (sigmoid * (qik - qjk) + reg * puk);
                Q(i,k) = Q(i,k) - lRate * (sigmoid * puk + reg * qik);
                Q(j,k) = Q(j,k) - lRate * (sigmoid * (-puk) + reg * qjk);
                loss = loss + reg * puk * puk + reg * qik * qik + reg * qjk * qjk;
            end
        end
    %calculate loss
    if abs(last_loss-loss)< 0.00001
        break;
    end
    last_loss=loss;
end
disp(steps)
eR = P*transpose(Q);
eR= transpose(eR);
%T1 = graythresh(eR);
%T= findThresh(eR);
T = graythresh(eR);
disp(T);
BW = im2bw(eR,T);
eR(BW==1)=2;
eR(BW==0)=1;
% eR(eR<=1+T)=1;
% eR(eR>1+T)=2;
dmiss(miss_idx)= eR(miss_idx);
eRc = zeros(size(data));
eRc(data==1 & eR==1) = 1;
eRc(data==2 & eR==2) = 2;
eRc(data==1 & eR==2) = 3;
eRc(data==2 & eR==1) = 4;

impute_res = zeros(size(data));
impute_res(data==1 & dmiss==1) = 1;
impute_res(data==2 & dmiss==2) = 2;
impute_res(data==1 & dmiss==2) = 3;
impute_res(data==2 & dmiss==1) = 4;



error_impute=size(find((data(miss_idx)-eR(miss_idx)) ~=0),1);
sprintf('total error: %d',size(eRc(eRc==4),1)+size(eRc(eRc==3),1))
sprintf('total one error: %d',size(eRc(eRc==3),1))
sprintf('total two error: %d',size(eRc(eRc==4),1))
sprintf('total imputes is: %d ', size(miss_idx,1))
sprintf('impute one error: %d',size(impute_res(impute_res==3),1))
sprintf('impute two error: %d',size(impute_res(impute_res==4),1))
sprintf('impute error is: %d ', error_impute)
figure,imagesc(data);title('original data');
figure,imagesc(dmiss);title('Result');
figure,imagesc(eRc);title('Error Encoded result');
figure,imagesc(impute_res);title('Error Encoded impute result');
end
