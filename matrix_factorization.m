function [P,Q] = matrix_factorization(R, P, Q, K,steps)
Errors = [];

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
     step
end
%save('P_replicate_slow.mat','P');
%save('Q_replicate_slow.mat','Q');
figure,plot(Errors);
end