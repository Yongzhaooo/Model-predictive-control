function [K,Pend]=DP(A,B,Q,R,Pf,N)
% A B is the state space matrix
% Q R Pf are the weight matrix
% N is the horizon steps
% P is the final cost
% k is the gain of U=Kx
P=cell(N+1,1);
P{N+1}=Pf;
for i = N:-1:1
P{i} = Q + A'*P{i+1}*A - A'*P{i+1}*B*inv(R + B'*P{i+1}*B)*B'*P{i+1}*A;
end
K= -inv(R+B'*P{2}*B)*B'*P{2}*A;
Pend=P{1};
end