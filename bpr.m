function eRc = bpr()
%[data,dmiss] = dataMiss(1,0.4);
data= load('cf_gt_data.mat');
data= data.data;
dmiss= load('cf_Miss_data.mat');
dmiss= dmiss.data_miss;

%make matrix 
users= size(dmiss,2);
items= size(dmiss,1);
figure,imagesc(dmiss);title(sum(sum(dmiss==0))/(items*users))
miss_idx = find(dmiss==0);
dmiss(dmiss==2)=-1;
tic;

K= 3;
P = load('P.mat');
P=P.P;
Q = load('Q.mat');
Q=transpose(Q.Q);
last_loss=0;
lRate= 0.001;
reg=0.001;
x = P* transpose(Q);           
% mat3d = zeros(items,items,users);
% loss= zeros(items,items,users);
% sigmoid= zeros(items,items,users);
% for u= 1:users
%     for i=1: items
%         for j=1 : items
%             mat3d(j,i,u)= sign(x(u,i)- x(u,j));
%             sigmoid(j,i,u)= 1/(1+exp(-1*mat3d(j,i,u)));
% 
%         end
%     end
% end
                
for steps =1:5000
    loss=0;
    for u= 1:users
        for i=1: items
            for j=1 : items
                if i ~=j
    %                 xui = predict(P,u,Q,i);
    %                 xuj = predict(P,u,Q,j);            
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
            end
        end
    end
    %calculate loss
    if abs(last_loss-loss)< 0.001
        break;
    end
    last_loss=loss;
end
eR = P*transpose(Q);
eR= transpose(eR);
T= findThresh(eR);
disp(T);
eR(eR<1+T)=1;
eR(eR>=1+T)=2;
dmiss(miss_idx)= eR(miss_idx);
eRc = zeros(size(data));
eRc(data==1 & dmiss==1) = 1;
eRc(data==2 & dmiss==2) = 2;
eRc(data==1 & dmiss==2) = 3;
eRc(data==2 & dmiss==1) = 4;

error_impute=size(find((data(miss_idx)-eR(miss_idx)) ~=0),1);
sprintf('total missing values: %d',size(miss_idx,1))
sprintf('impute error is: %d ', error_impute)
figure,imagesc(data);title('original data');
figure,imagesc(dmiss);title('Result');
figure,imagesc(eRc);title('Error Encoded result');
time=toc;
disp(time);
end

% function x= predict(P,u,Q,i)
% x = P(u,:) * transpose(Q(i,:));
% end
