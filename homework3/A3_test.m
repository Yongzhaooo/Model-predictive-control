clear;
clc;
% Define A and b
A = [0.4889 0.2939; 1.0347 -0.7873; 0.7269 0.8884; -0.3034 -1.1471];
b = [-1.0689; -0.8095; -2.9443; 1.4384];

% Reshape A into a matrix of size 2n by n
F = [A, -ones(4,1); -A, -ones(4,1)];
g = [b; -b];
c = [0 0 1];

% Solve the linear program
x = linprog(c, F, g)

% Extract the solution for x
x_sol = x(1:2);

% Print the solution
disp(x_sol);
y=A*x_sol-b


%% 4 Quadratic programming
% expand x to [x1 x2 u0 u1]
H = eye(4);
f = zeros(4,1);

% xk+1=0.4xk+uk can be translate into 2 euqations:
% x1 - u0 = 1.5
% x2 - 0.4x1 -u1 = 0

Aeq = [1 0 -1 0;
       -0.4 1 0 -1];
beq = [1.5;0];

%  dispate the inequality to single side for the standard form
Aleq= [eye(4);
       -eye(4)];
bleq= [5 0.5 2 2   -2.5 0.5 2 2]';

% lb=[2.5 -0.5 -2 -2]';
% ub=[5 0.5 2 2]'

options = optimoptions('quadprog','Display','iter');
[x,fval,exitflag,output,lambda] = quadprog(H,f,Aleq,bleq,Aeq,beq,[],[],[],options);
% [x,fval,exitflag,output,lambda] = quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);

% Print the solution and other relevant information
disp(x);
disp(fval);

%verify the constarins
mu = lambda.ineqlin
if all(mu >= 0)
disp('mu >= 0 is satisfied');
else
    disp('the constraint for mu is not satisfied');
end

gx=Aleq*x-bleq
if all(gx <= 0)
    disp('gx <= 0 is satisfied');
else
    disp('the constraint for gx is not satisfied');
end

hx=Aeq*x-beq
if all( abs(hx) <= 1e-5)
    disp('hx == 0 is satisfied');
else
    disp('the constraint for hx is not satisfied');
end

muigi=mu.*gx
if all(abs(muigi) <= 1e-5)
    disp('muigi == 0 is satisfied');
else
    disp('the constraint for muigi is not satisfied');
end

kkt_residual=x + Aleq.'*lambda.ineqlin + Aeq'*lambda.eqlin

if norm(kkt_residual) <= 1e-5
    disp('KKT residual is satisfied');
else
    disp('KKT residual is not satisfied');
end


%% Remove the x1 lower bound
Aleq1 = Aleq([1:4],:);
bleq1 = bleq(1:4);

options = optimoptions('quadprog','Display','iter');
[x,fval,exitflag,output,lambda] = quadprog(H,f,Aleq1,bleq1,Aeq,beq,[],[],[],options)

x
fval
lambda.ineqlin

%% Remove the x1 upper bound

Aleq2 = Aleq([5:8],:);
bleq2 = bleq(5:8);

options = optimoptions('quadprog','Display','iter');
[x,fval,exitflag,output,lambda] = quadprog(H,f,Aleq2,bleq2,Aeq,beq,[],[],[],options)

x
fval
lambda.ineqlin
