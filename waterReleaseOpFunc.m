function [f,FE] = waterReleaseOpFunc(Q,Demand,FE)
% objective function (minimization)
%f=sum((Q-Demand).^2);

[M,N] = size(Q);
f = 0;
for i = 1: M
    for j = 1:N
       f = f+ ((Q(i,j)- Demand(i,j))^2) ;
    end
end
FE = FE + 1;
end

