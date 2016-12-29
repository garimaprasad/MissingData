noise=[0.05,0.2,0.4,0.6,0.8];
%noise=[0.2,0.3];
r=size(noise,2);
temp=[1,3,5];
replicate=5;
I_1_err=zeros(r,replicate);
I_err=zeros(r,replicate);
T_err=zeros(r,replicate);
I_2_err=zeros(r,replicate);
M_one=zeros(r,replicate);
M_two=zeros(r,replicate);
N_per=zeros(r,replicate);
R_1=zeros(r,replicate);
R_2=zeros(r,replicate);
R_1_err=zeros(r,replicate);
R_2_err=zeros(r,replicate);
dorep=false;
%ground truth data
data= load('data.mat');
data= data.data;
x = data;
K=3;
%replicate data
for i = 1:replicate
    if dorep
        data = cat(2,x,data);
    end
    dorep=true;
    for n=1:r
        [rem_1,rem_2,miss_one,miss_two,n_per,imp_1_err,imp_2_err,rem_1_err,rem_2_err,error_impute]=test_cf(data,noise(1,n),K);
        I_1_err(n,i)=imp_1_err;
        I_2_err(n,i)=imp_2_err;
        M_one(n,i)=miss_one;
        M_two(n,i)=miss_two;
        N_per(n,i)=n_per;
        R_1(n,i)=rem_1;
        R_2(n,i)=rem_2;
        R_1_err(n,i)=rem_1_err;
        R_2_err(n,i)=rem_2_err;
        T_err(n,i)=(imp_1_err+imp_2_err+rem_1_err+rem_2_err)*100/(size(data,1)*size(data,2));
        I_err(n,i)= (error_impute/(miss_one+ miss_two))*100;
    end
    
end

I_1_err= (I_1_err./M_one)*100;
I_2_err= (I_2_err./M_two)*100;
R_1_err=(R_1_err./R_1)*100;
R_2_err=(R_2_err./R_2)*100;

figure;
for j = temp 
    %axis([0 50 0 30])
    plot(N_per(:,j),I_1_err(:,j),'o-', 'DisplayName',num2str(j))
    title('one error on the imputes')
    xlabel('noise percent') % x-axis label
    ylabel('error percent') % y-axis label
    hold on;
end
legend('show')
hold off;
figure;
for j = temp 
    %axis([0 50 0 30])
    plot(N_per(:,j),I_2_err(:,j),'o-' ,'DisplayName',num2str(j))
    title('two error on the imputes')
    xlabel('noise percent') % x-axis label
    ylabel('error percent') % y-axis label
    hold on;
end
legend('show')
hold off;
figure;
for j = temp
    %axis([0 50 0 30])
    plot(N_per(:,j),R_1_err(:,j),'o-','DisplayName',num2str(j))
    title('one error on the remaining values')
    xlabel('noise percent') % x-axis label
    ylabel('error percent') % y-axis label  
    hold on;
end
legend('show')
hold off;
figure;
for j = temp
    %axis([0 50 0 30])
    plot(N_per(:,j),R_2_err(:,j),'o-','DisplayName',num2str(j))
    title('two error on the remaining values')
    xlabel('noise percent') % x-axis label
    ylabel('error percent') % y-axis label  
    hold on;
end
legend('show')
hold off;
figure;
for j = temp
    %axis([0 50 0 30])
    plot(N_per(:,j),I_err(:,j),'o-','DisplayName',num2str(j))
    title('impute error')
    xlabel('noise percent') % x-axis label
    ylabel('error percent') % y-axis label
    
    hold on;
end
legend('show')
hold off;
figure;
for j = temp 
    
    %axis([0 50 0 30])
    plot(N_per(:,j),T_err(:,j),'o-','DisplayName',num2str(j))
    title('Total error')
    xlabel('noise percent') % x-axis label
    ylabel('error percent') % y-axis label
    
    hold on;
end
legend('show')
hold off;
plot(I_2_err);
