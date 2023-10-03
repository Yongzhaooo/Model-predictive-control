function [Z,VN]=CRHC(A,B,N,Q,R,Pf,F,G,h,x0,n)
%% This function solves a quadratic program (QP) to compute the optimal control 
% of a constrained finite-horizon LQ problem, which is in uncondensed form,
% i.e. both x and u are optimization variables.
%% A and B are the system matrices when x(k+1)=Ax(k)+Bu(k)
%% Q, R, and Pf are the gains in the cost function
%% N is the length of the horizon
%% Z is the vector of optimal variables and VN is the cost function 
%% F, G, h are constraint matrices
%% The constraints have the form: Fx + Gu <= h  
%% x0 is the initial condition


% Objective function
Qbar = blkdiag(kron(eye(N-1),Q),Pf);
Rbar = kron(eye(N),R);
H    = blkdiag(Qbar,Rbar);
f    = [];

% Equality constraints
I    = eye(n);
Aeq1 = kron(eye(N),I)+kron(diag(ones(N-1,1),-1),-A);
Aeq2 = kron(eye(N),-B);
Aeq  = [Aeq1 Aeq2];
beq  = [A*x0;zeros(n*(N-1),1)];

% Inequality constraints
Ain  = [F G];
bin  = h;

% Solve QP
[Z,VN,exitflag,output,lambda] = quadprog(2*H,f,Ain,bin,Aeq,beq);

% Check if the problem has been solved to optimality
if exitflag ~= 1
    disp(['The exitflag is ' num2str(exitflag) ', which means that something is wrong.'])
    pause
end


% Note: bounds on the optimization variables can be passed directly to 
% quadprog (lb,ub) instead of actual inequality constraints.

% Note: one could also use the condensed version. In that case, the only 
% optimization variables of the QP would be the control inputs. 
