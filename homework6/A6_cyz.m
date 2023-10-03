%part1 Lyapunov functions

clc
clear
%define Matrices
A = [1.0041 0.0100 0 0;
    0.8281 1.0041 0 -0.0093;
    0.0002 0.0000 1 0.0098;
    0.0491 0.0002 0 0.9629];
B = [0.0007;
    0.1398;
    0.0028;
    0.5605];

%part(a)
S=eye(4);
Q=S-A'*S*A;
disp('-----Q matrix------')
disp(Q);

%validate
Q'-Q

% eignQ=eig(Q)
% if all((eignQ)>0)
%     disp('Q has all positve eign values')
% else
%     disp('Q has negative eigen value')
% end
% 
% eigenA = eig(A)
% if all(abs(eigenA) < 1)
%     disp("All eigenvalues of A is smaller than 1")
% else
%     disp("Not all eigenvalues of A are smaller than 1")
% end

eigQ=checkQ(Q);
eigA=checkA(A);

%partb
K = [114.3879 12.7189 -1.2779 -1.5952];
%Feedback Ak
Ak = A-B*K;
Qb = S-Ak'*S*Ak
eigQak=checkA(Ak)
eigQb=checkQ(Qb)

%partc
Sc = dlyap(Ak',Q)
issymmetric(Sc)
% checkQ(Qc)
checkQ(Sc)

disp('-------part 2----------')
%2. Stability with receding horizon control
clear
A = [1.0041 0.0100 0 0;
    0.8281 1.0041 0 -0.0093;
    0.0002 0.0000 1 0.0098;
    0.0491 0.0002 0 0.9629];
B = [0.0007;
    0.1398;
    0.0028;
    0.5605];
R = 1;
Q = eye(4);
%2a
%2b
Pf = Q;
for N = 1:100
    [P0,K] = DP(A,B,Pf,Q,R,N);
    Ak = A-B*K;
    if all(abs(eig(Ak))<1)
        disp("The shortest N is: ")
        disp(N);
        disp("The K is")
        disp(K);
        break;
    end
end
P0
%2c
Qc=P0-(A+B*K)'*P0*(A+B*K)
% issymmetric(Qc)
[~,flag] = chol(Qc)
checkQ(Qc)
%d
[Pf,~,K] = dare(A,B,Q,R)
N = 1;
% [P0,K] = DP(A,B,Pf,Q,R,N); % same K as above
Ak = A-B*K;
Qc = Ak'*Pf*Ak-Pf+Q+(K'*B'*R*B*K)
[~,flag] = chol(Qc)
checkQ(Qc)

function eignQ = checkQ(Q)

eignQ = eig(Q);  % 计算Q的特征值

if all((eignQ) > 0)
    disp(sprintf("All eig of %s is positive.", inputname(1)));  % 如果所有特征值都为正，输出系统稳定
else
    disp(sprintf("Not all eigenvalues of %s are positive.", inputname(1)));  % 如果有负的特征值，输出系统不稳定
end

end


function eigenA = checkA(A)

eigenA = eig(A);  % 计算A的特征值

if all(abs(eigenA) < 1)
   disp(sprintf("All eig of %s is smaller than 1.", inputname(1)));  % 如果所有特征值的绝对值都小于1，输出系统稳定
else
    disp(sprintf("Not all eigenvalues of %s are smaller than 1",inputname(1)))  % 如果有特征值的绝对值大于等于1，输出系统不稳定
end

end

function [P0,K] = DP(A,B,Pf,Q,R,N)
    P = cell(1,N+1);
    P{N+1} = Pf;
    for i = N+1:-1:2
        P{i-1} = Q + A'*P{i}*A - A'*P{i}*B*inv(R + B'*P{i}*B)*B'*P{i}*A;
    end
    K = inv(R+B'*P{2}*B)*B'*P{2}*A;
    P0 = P{1};
end

