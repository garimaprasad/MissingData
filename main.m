%find best k
synthData(1000,10);
error= zeros(1000);
time= zeros(1000);

    for k = 1:1000
        tic
        error(k)= cf(50,0.2,k);
        time(k)= toc;
        disp(time(r,k))
    end
figure,plot(error);
[minvalue, minidx] = min(error);
k=minidx;

