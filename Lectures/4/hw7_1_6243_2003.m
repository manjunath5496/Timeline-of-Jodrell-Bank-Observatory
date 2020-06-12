function hw7_1_6243_2003(p,r)
% function hw7_1_6243_2003(p,r)
%
% searches for unstable unimodal cycles of half-period 1 for dx/dt=Ax+Bsign(Cx),
% where p is the characteristic polynomial of  A 

[A,B]=ssdata(tf(1,p));
eA=expm(A);
n=size(A,1);
I=eye(n);
x0=inv(I+eA)*(I-eA)*inv(A)*B;
C=[0 1 r]/[x0 B A*x0];

N=300;
t=linspace(0,1,N);
y=lsim(ss(A,B,C,0),ones(N,1),t,x0);
y=squeeze(y);
abs(eig(eA-(B-A*x0)*C*eA/(C*(B-A*x0))))
close(gcf); plot(t,y); grid
    