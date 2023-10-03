%% Main script of PSS1 in SSY281
clear;
clc;
close all;

A =[
    1.0041, 0.0100, 0, 0;
    0.8281, 1.0041, 0, -0.0093;
    0.0002, 0.0000, 1, 0.0098;
    0.0491, 0.0002, 0, 0.9629];

B =[
    0.0007;
    0.1398;
    0.0028;
    0.5605];

C=[1,0,0,0;0,0,1,0];

Q=eye(4);

Pf=kron(eye(4),10)

% R=1;

x0=[pi/38,0,0,0]'; % initial state

n  = length(A);
m  = size(B,2);

% N  = 10;    % set horizon length

R4=[1,1,0.1,0.1];
N4=[40,80,40,80];

Kc=cell(4,1);
x_vec  = cell(4,1); % store state trajectory
u_vec  = cell(4,1);      % store control input trajectory

for i=1:4
    R=R4(i);
    N=N4(i);
    % Assemble constraints in matrix form as: Fx + Gu <= h
    % Here, assume the box constraints: abs(x2(k)) <= x2_max and abs(u(k)) <= u_max for all k
    x2_max = 1;
    u_max  = 8;
    F      = kron([eye(N); zeros(N)], [0 1 0 0;0 -1 0 0;0 0 0 0;0 0 0 0]);
    G      = kron([zeros(N); eye(N)], [1; -1; 1 ; -1]);
    h      = [x2_max*ones(n*N,1); u_max*ones(n*N,1)];
    % Simulation parameters
    tf     = 200;               % assume a sampling interval of 1
    x      = x0;               % initial condition

    % Carry out simulations by applying the receding horizon principle
    for iter = 1:tf
%         % store state and control input for each iteration
%         x_vec{i}  = [x0 zeros(n,tf)]; % store state trajectory
%         u_vec{i}  = zeros(m,tf);      % store control input trajectory
        [Z,~] = CRHC(A,B,N,Q,R,Pf,F,G,h,x,n);   % constrained RHC
        % Get next command and next state
        x     = Z(1:n);
        u     = Z(n*N+1);
        x_vec{i}(:,iter+1) = x;
        u_vec{i}(:,iter)   = u;
    end
%     x_vec
%     u_vec
end

figure(1)
hold on
grid on
% Subplot 1
subplot(4,1,1)
ax1 = gca; % Save the current axis handle
xlabel('time-step k'), ylabel('x1(k)')
step1 = 2:1:101;
legend_strings1 = cell(1,4);
for i = 1:4
    plot(step1,x_vec{i}(1,2:101),'LineWidth',2);
    legend_strings1{i} = sprintf('x(1) R=%.2f N=%d',R4(i), N4(i));
    hold on
end
legend(ax1,legend_strings1, 'Location', 'best');


% Subplot 2
subplot(4,1,2)
ax2 = gca;
xlabel 'time-step k', ylabel 'x2(k)'
step2=2:1:101;
legend_strings2 = cell(1,4);
for i =1 :4
    plot(step2,x_vec{i}(2,2:101),'LineWidth',2);
    legend_strings2{i} = sprintf('x(2) R=%.2f N=%d',R4(i), N4(i));
    hold on
end
legend(ax2,legend_strings2, 'Location', 'best');


% Subplot 3
subplot(4,1,3)
ax3 = gca;
xlabel 'time-step k', ylabel 'x3(k)'
step3=2:1:201;
legend_strings = cell(1,4);
for i =1 :4
    plot(step3,x_vec{i}(3,2:201),'LineWidth',2);
    legend_strings{i} = sprintf('x(3) R=%.2f N=%d',R4(i), N4(i));
    hold on
end
legend(ax3,legend_strings, 'Location', 'best');

% Subplot 4
subplot(4,1,4)
ax4 = gca;
xlabel 'time-step k', ylabel 'u(k)'
step4=2:1:31;
legend_strings = cell(1,4);
for i =1 :4
    plot(step4,u_vec{i}(2:31),'LineWidth',2);
    legend_strings{i} = sprintf('u(k) R=%.2f N=%d',R4(i), N4(i));
    hold on
end
legend(ax4, legend_strings, 'Location', 'best');






