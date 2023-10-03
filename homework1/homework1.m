clear;
clc;
close all;

A =[
1.0041, 0.0100, 0, 0;
0.8281, 1.0041, 0, -0.0093;
0.0002, 0.0000, 1, 0.0098;
0.0491, 0.0002, 0, 0.9629]

B =[
0.0007;
0.1398;
0.0028;
0.5605]

C=[1,0,0,0;0,0,1,0]

Q=eye(4)

Pf=kron(eye(4),10)

R=1

%for CHRC
n  = length(A); 
m  = size(B,2);
x0=[pi/38,0,0,0]'; % initial state

%(a)1 DP method
% set a flag
stable = false;
Na = 0;
while ~stable
    Na = Na+1;
    [Kend,Pend] = DP(A,B,Q,R,Pf,Na);
    if all(abs(eig(A+B*Kend))<=1)
        stable = true;
    end
end
Na
eiga=eig(A+B*Kend)
Pend
Kend

%(b)P_inf
[P,F,G] = idare(A,B,Q,R)
[Kcon,Pcon,timecon]=DP_con(A,B,Q,R,Pf,1000)

%(c) 
stable = false;
Nc = 0;
while ~stable
    Nc = Nc+1;
    [Kend1,Pend1] = DP(A,B,Q,R,Pcon,Nc);
    if all(abs(eig(A+B*Kend))<=1)
        stable = true;
    end
end

Nc
eigc=eig(A+B*Kend1)
Pend1
Kend1

% 应该是因为第二步已经算出来了静态解，所以直接就收敛了

%3  batch solution
stable = false;
N5 = 0;
while ~stable
    N5 = N5+1;
    [K0,Pend] = mBatch(A,B,Q,R,Pf,N5);
    if all(abs(eig(A+B*K0))<=1)
        stable = true;
    end
end

Pf
N5
eig(A+B*K0)
Pend
K0


%4. Receding horizon control
R4=[1,1,0.1,0.1];
N4=[40,80,40,80];
Kend=cell(4,1);
%using mBatch to solve this
for i = 1:4
    [Kend{i},Pend]=mBatch(A,B,Q,R4(i),Pf,N4(i));
end

Kend
T=201;
x = cell(4,1);
u = cell(4,1);
for i = 1:4
    x{i} = zeros(T, 4);
    u{i} = zeros(T, 1);
    x{i}(1, :) = x0;
    for j = 2:T
        x{i}(j, :) = (A + B * Kend{i}) * x{i}(j-1, :)';
        u{i}(j) = Kend{i}*x{i}(j-1, :)';
    end
end

% plot part
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
    plot(step1,x{i}(2:101,1),'LineWidth',2);
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
    plot(step2,x{i}(2:101,2),'LineWidth',2);
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
    plot(step3,x{i}(2:201,3),'LineWidth',2);
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
    plot(step4,u{i}(2:31),'LineWidth',2);
    legend_strings{i} = sprintf('u(k) R=%.2f N=%d',R4(i), N4(i));
    hold on
end
legend(ax4, legend_strings, 'Location', 'best');

%fuck this legend work, it waste me about an hour to do
