clc;
clear;

%%1
a = -0.9421; b = 82.7231; c = 14.2306; p = -3.7808; q = 4.9952; r = 57.1120;
h= 0.1;
A=[0,1,0,0;
    b,0,0,a;
    0,0,0,1;
    q,0,0,p]
B=[0;c;0;r]
C=[1,0,0,0;0,0,1,0]
D=0

disH=expm([A,B;zeros(1,4),0]*h)
Adh=disH(1:4,1:4)
Bdh=disH(1:4,5)

%verify
sysc=ss(A,B,C,D);
sysd=c2d(sysc,h)

%%2 
t=0.8*h;
disHt=expm([A,B;zeros(1,4),0]*(h-t))
Adht=disHt(1:4,1:4)
Bdht=disHt(1:4,5)

% eig(A)
eig(Adh)
eig(Adht)
