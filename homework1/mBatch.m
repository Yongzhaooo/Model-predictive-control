function [Kend,Pend]=mBatch(A,B,Q,R,Pf,N)
% A B is the state space matrix
% Q R Pf are the weight matrix
% P is the final cost
% k is the gain of U=Kx
% N is the horizon steps
gamma = kron(eye(N),B);
omega = A;
for i=1:N-1
    gamma = gamma + kron(diag(ones(N-i,1),-i),A^i*B);
    omega = [omega; A^(i+1)];
end

Qb = blkdiag( kron(eye(N-1),Q), Pf );
Rb = kron(eye(N),R);

Pend = (Q + omega'*Qb*omega - omega'*Qb*gamma*inv(gamma'*Qb*gamma + Rb)*gamma'*Qb*omega);

K = -inv(gamma'*Qb*gamma + Rb)*gamma'*Qb*omega;
Kend = K(1,:);
end
