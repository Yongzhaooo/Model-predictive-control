function [K,Pend,time]=DP_con(A,B,Q,R,Pf,N)
% A B is the state space matrix
% Q R Pf are the weight matrix
% P is the final cost
% k is the gain of U=Kx
% N is the horizon steps
P=cell(N+1,1);
P{N+1}=Pf;
for i = N:-1:1
    P{i} = Q + A'*P{i+1}*A - A'*P{i+1}*B*inv(R + B'*P{i+1}*B)*B'*P{i+1}*A;
    if norm(P{i+1}-P{i}) <= 1e-1
        break;
    end
end
if i==1
    disp('without con')
else
    K= -inv(R+B'*P{i+1}*B)*B'*P{i+1}*A;
end
Pend=P{i};
time=N-i+1;
end