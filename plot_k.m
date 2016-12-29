
k = [2,6,10,16,20,30];
%r=[1, 2, 5,10];
n = ['r','l','s', 'u','b' ];
r=ones(length(n));
total=[];
impute=[];
for i= 1:5
    fname= sprintf('/Users/Garima/Documents/8980DM/project/dataset/Garima2/%s1h99v99_1.mat',n(i));
    load(fname);
    mapStack= mapStack(:,1:r(i)*100);

    static_inds = sum(mapStack==1,2)==r(i)*100 | sum(mapStack==2,2)==r(i)*100;
    %static_inds = sum(mapStack==1,2)==100 | sum(mapStack==2,2)==100;
    R = mapStack(static_inds==0,:);
    GT= GT(:,1:r(i)*100);
    data= GT(static_inds==0,:);
    miss_idx = find(R==0);
    remain_idx= find(R>0);
    
    path= sprintf('/Users/Garima/Documents/8980DM/project/final_data_results/%sP_16.mat',n(i));
    P1= load(path);
    P1= P1.P;
    path= sprintf('/Users/Garima/Documents/8980DM/project/final_data_results/%sQ_16.mat',n(i));
    Q1= load(path);
    Q1= Q1.Q;
    
    eR= P1*Q1;

    T= findThresh(P1,Q1,remain_idx,R);
    eR(eR<=T)=1;
    eR(eR>T)=2;
    R(miss_idx)= eR(miss_idx);

    eRc = zeros(size(data));
    eRc(data==1 & eR==1) = 1;
    eRc(data==2 & eR==2) = 2;
    eRc(data==1 & eR==2) = 3;
    eRc(data==2 & eR==1) = 4;

    impute_res = zeros(size(data));
    impute_res(data==1 & R==1) = 1;
    impute_res(data==2 & R==2) = 2;
    impute_res(data==1 & R==2) = 3;
    impute_res(data==2 & R==1) = 4;

    total_err= (size(eRc(eRc==3),1) +size(eRc(eRc==4),1))*100/numel(data)
    imp_err= (size(impute_res(impute_res==3),1)+size(impute_res(impute_res==4),1))*100/numel(miss_idx)
    total= [total total_err];
    impute= [impute imp_err];
end


k= [2,6,8,12,16,20,30];
total=[ 6.2668,3.2518,3.3104,2.9862,3.6578,3.4183,3.2601];
impute=[11.8508,10.5660,9.6884,11.1172,15.9254,16.9133,20.5130];
figure;
plot(k1,t1,'o-')
title('Total Error for Spatio-Temporal holes and noise')
xlabel('K') % x-axis label
ylabel('Error Percent') % y-axis label

figure;
plot(k1,i1,'o-')
title('Error for holes for Spatio-Temporal holes and noise')
xlabel('K') % x-axis label
ylabel('Error Percent') % y-axis label
