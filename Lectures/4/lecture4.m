%
% B747 longitudinal dynamics
% AA271a
% B747 at Mach 0.8, 40,000ft, level-flight
% From Etkin and Reid page 166
%
if 1
	% metric   
    Xu=-1.982e3;Xw=4.025e3;
	Zu=-2.595e4;Zw=-9.030e4;Zq=-4.524e5;Zwd=1.909e3;
	Mu=1.593e4;Mw=-1.563e5;Mq=-1.521e7;Mwd=-1.702e4;

	g=9.81;theta0=0;S=511;cbar=8.324;
	U0=235.9;Iyy=.449e8;m=2.83176e6/g;cbar=8.324;rho=0.3045;
	Xdp=.3*m*g;Zdp=0;Mdp=0;
	Xde=-3.818e-6*(1/2*rho*U0^2*S);Zde=-0.3648*(1/2*rho*U0^2*S);
	Mde=-1.444*(1/2*rho*U0^2*S*cbar);
   
else
	%english   
	Xu=-1.358e2;Xw=2.758e2;cbar=27.31;
	Zu=-1.778e3;Zw=-6.188e3;Zq=-1.017e5;Zwd=1.308e2;
	Mu=3.581e3;Mw=-3.515e4;Mq=-1.122e7;Mwd=-3.826e3;

	g=32.2;U0=774;Iyy=.331e8;m=636636/g;cbar=27.31;rho=0.0005909;
	theta0=0;S=5500;
	Xdp=.3*m*g;Zdp=0;Mdp=0;
	Xde=-3.818e-6*(1/2*rho*U0^2*S);Zde=-0.3648*(1/2*rho*U0^2*S);
	Mde=-1.444*(1/2*rho*U0^2*S*cbar);
   
end

A=[Xu/m Xw/m 0 -g*cos(theta0);[Zu Zw Zq+m*U0 -m*g*sin(theta0)]/(m-Zwd);
   [Mu+Zu*Mwd/(m-Zwd) Mw+Zw*Mwd/(m-Zwd) Mq+(Zq+m*U0)*Mwd/(m-Zwd) ...
         -m*g*sin(theta0)*Mwd/(m-Zwd)]/Iyy;
   [ 0 0 1 0]];

B=[Xde/m Xdp/m;Zde/(m-Zwd) Zdp/(m-Zwd);(Mde+Zde*Mwd/(m-Zwd))/Iyy ...
      (Mdp+Zdp*Mwd/(m-Zwd))/Iyy;0 0];
C=[eye(3,4)];

wsp=sqrt(-U0*Mw/Iyy)
zetasp=-Mq/2*sqrt(-1/U0/Mw/Iyy)

wsp=sqrt(Zw*Mq/m/Iyy-U0*Mw/Iyy)
zetasp=-(Zw/m+Mq/Iyy+Mwd/Iyy*U0)/2/wsp

[V,ev]=eig(A);ev=diag(ev);D=V;
%rifd(ev)

% nondimensionalization of the eigenvetors
V(1,:)=V(1,:)/U0;
V(2,:)=V(2,:)/U0;
V(3,:)=V(3,:)/(2*U0/cbar);
for i=1:4
   V(:,i)=V(:,i)/V(4,i);
end

figure(3)
subplot(221)
polar(phase(V(:,1)),abs(V(:,1)))
axis('square')
xlabel(num2str(abs(ev(1))))
subplot(222)
polar(phase(V(:,3)),abs(V(:,3)))
axis('square')
xlabel(num2str(abs(ev(3))))
%
% initial condition response using the original EVectors
% as the the IC
%
sys=ss(A,zeros(size(A,1),1),eye(3,4),zeros(3,1));
[y1,t1,x1]=initial(sys,real(D(:,1)));
[y3,t3,x3]=initial(sys,real(D(:,3)));

subplot(223)
plot(t1,y1(:,1),'-',t1,y1(:,2),'--',t1,y1(:,3),'-.')
ylabel('Perturbation States u,w,q');xlabel('time (sec)');axis([0 15 -1 .5])
subplot(224)
plot(t3,y3(:,1),'-',t3,y3(:,2),'--',t3,y3(:,3),'-.')
ylabel('Perturbation States u,w,q');xlabel('time (sec)');axis([0 600 -1 1])

%return
%      real       imaginary     frequency      damping
% 
%  -3.2889e-003   -6.7202e-002    6.7282e-002    4.8882e-002 
%  -3.7168e-001    8.8692e-001    9.6166e-001    3.8650e-001 
%V 
%   0.0156 + 0.0244i    -0.0254 - 0.6165i
%   1.0202 + 0.3553i     0.0045 - 0.0356i
%  -0.0066 + 0.0156i    -0.0001 - 0.0012i
%   1.0000               1.0000          

%
% step response to elevator - scaled to del_e=1 deg input
%
C=[eye(4,4);0 -1/U0 0 1];
sys_sc=ss(A,[B(:,1)/180*pi],C,zeros(5,1));
[Y,T]=step(sys_sc);
[Y2,T2]=step(sys_sc,[0:.1:40]);
%
% step response to thrust - scaled to del_p=1/6 
%
sys_sc=ss(A,[B(:,2)/6],C,zeros(5,1));
[Yt,Tt]=step(sys_sc);
[Yt2,Tt2]=step(sys_sc,[0:.1:40]);
%
%final value theorem:
Yfinal=-C*inv(A)*[B(:,1)/180*pi B(:,2)/6];
Yfinal(2,:)=Yfinal(2,:)/U0 % switch the w to w/U0=alpha
%
% freq response
%
freq=logspace(-2,1,400);
sys=ss(A,B,C,zeros(5,2));
G=freqresp(sys,freq);

clear Gude Gwde Gtde Gudp Gwdp Gtdp Ggde Ggdp
Gude(:,1)=G(1,1,:);
Gwde(:,1)=G(2,1,:);
Gtde(:,1)=G(4,1,:);
Ggde(:,1)=G(5,1,:);
Gudp(:,1)=G(1,2,:);
Gwdp(:,1)=G(2,2,:);
Gtdp(:,1)=G(4,2,:);
Ggdp(:,1)=G(5,2,:);

figure(1);clf
orient('landscape')
subplot(231)
loglog(freq,abs(Gude),min(abs(ev))*[1 1],[.1 1e4],'--',max(abs(ev))*[1 1],[.1 1e4],'--')
axis([.01 5 1 1e4]) 
ylabel('|G_u_d_e|')
xlabel('Freq (rad/sec)')
subplot(232)
loglog(freq,abs(Gwde)/U0,min(abs(ev))*[1 1],[.1 1e4],'--',max(abs(ev))*[1 1],[.1 1e4],'--')
axis([.01 5 .1 1e4]) 
ylabel('|G_{\alpha_d_e}|')
xlabel('Freq (rad/sec)')
title('Transfer function from elevator to flight variables')
subplot(233)
loglog(freq,abs(Ggde),min(abs(ev))*[1 1],[.1 1e4],'--',max(abs(ev))*[1 1],[.1 1e4],'--')
ylabel('|G_{\gamma_d_e}|')
xlabel('Freq (rad/sec)')
axis([.01 5 .1 1e4]) 
subplot(234)
semilogx(freq,180/pi*phase(Gude.'),min(abs(ev))*[1 1],[-360 0],'--',max(abs(ev))*[1 1],[-360 0],'--')
ylabel('arg G_u_d_e')
xlabel('Freq (rad/sec)')
axis([.01 5 -360 0]) 
subplot(235)
semilogx(freq,180/pi*phase(Gwde.')-360,min(abs(ev))*[1 1],[-360 0],'--',max(abs(ev))*[1 1],[-360 0],'--')
ylabel('arg G_{\alpha_d_e}')
xlabel('Freq (rad/sec)')
axis([.01 5 -360 0]) 
subplot(236)
semilogx(freq,180/pi*phase(Ggde.'),min(abs(ev))*[1 1],[-360 0],'--',max(abs(ev))*[1 1],[-360 0],'--')
ylabel('arg G_{\gamma_d_e}')
xlabel('Freq (rad/sec)')
axis([.01 5 -360 0]) 
print -depsc ac_fig1

figure(2);clf
orient('landscape')
subplot(231)
loglog(freq,abs(Gudp),min(abs(ev))*[1 1],[.01 1e4],'--',max(abs(ev))*[1 1],[.01 1e4],'--')
axis([.01 5 .01 1e3]) 
ylabel('|G_u_d_p|')
xlabel('Freq (rad/sec)')
subplot(232)
loglog(freq,abs(Gwdp)/U0,min(abs(ev))*[1 1],[.01 1e4],'--',max(abs(ev))*[1 1],[.01 1e4],'--')
axis([.01 5 .01 1e3]) 
ylabel('|G_{\alpha_d_p}|')
xlabel('Freq (rad/sec)')
title('Transfer function from thrust to flight variables')
subplot(233)
loglog(freq,abs(Ggdp),min(abs(ev))*[1 1],[.01 1e4],'--',max(abs(ev))*[1 1],[.01 1e4],'--')
ylabel('|G_{\gamma_d_p}|')
xlabel('Freq (rad/sec)')
axis([.01 5 .01 1e3]) 
subplot(234)
semilogx(freq,180/pi*phase(Gudp.'),min(abs(ev))*[1 1],[-360 360],'--',max(abs(ev))*[1 1],[-360 360],'--')
ylabel('arg G_u_d_p')
xlabel('Freq (rad/sec)')
axis([.01 5 -180 180]) 
subplot(235)
semilogx(freq,180/pi*phase(Gwdp.'),min(abs(ev))*[1 1],[-360 360],'--',max(abs(ev))*[1 1],[-360 360],'--')
ylabel('arg G_{\alpha_d_p}')
xlabel('Freq (rad/sec)')
axis([.01 5 -180 180]) 
subplot(236)
semilogx(freq,180/pi*phase(Ggdp.'),min(abs(ev))*[1 1],[-360 360],'--',max(abs(ev))*[1 1],[-360 360],'--')
ylabel('arg G_{\gamma_d_p}')
xlabel('Freq (rad/sec)')
axis([.01 5 -360 0]) 
print -depsc ac_fig2
%
% plot time stuff
%
figure(3);clf;orient('landscape')
LL1=1:400;LL2=1:350;
subplot(231)
plot(T(LL2),Y(LL2,1,1),'LineWidth',2)
ylabel('u');xlabel('time')
hold on;plot([400 600],[1 1]*Yfinal(1,1),'--');hold off
axis([0 600 0 30])
subplot(232)
plot(T(LL2),Y(LL2,2,1)/U0,'LineWidth',2)
ylabel('\alpha (rad)');xlabel('time')
hold on;plot([400 600],[1 1]*Yfinal(2,1),'--');hold off
axis([0 600 -.03 0])
title('Step response to 1 deg elevator perturbation')
subplot(233)
plot(T(LL2),Y(LL2,5,1),'LineWidth',2)
ylabel('\gamma');xlabel('time')
hold on;plot([400 600],[1 1]*Yfinal(5,1),'--');hold off
axis([0 600 -.1 .1])
subplot(234)
plot(T2(LL1),Y2(LL1,1,1),'LineWidth',2)
ylabel('u');xlabel('time')
subplot(235)
plot(T2(LL1),Y2(LL1,2,1)/U0,'LineWidth',2)
ylabel('\alpha (rad)');xlabel('time')
subplot(236)
plot(T2(LL1),Y2(LL1,5,1),'LineWidth',2)
ylabel('\gamma');xlabel('time')
print -depsc ac_fig3
%
%
%
figure(6);clf;orient('landscape')
LL1=1:400;LL2=1:350;
subplot(231)
plot(Tt(LL2),Yt(LL2,1,1),'LineWidth',2)
ylabel('u');xlabel('time')
hold on;plot([400 600],[1 1]*Yfinal(1,2),'--');hold off
axis([0 600 -15 15])
subplot(232)
plot(Tt(LL2),Yt(LL2,2,1)/U0,'LineWidth',2)
ylabel('\alpha (rad)');xlabel('time')
hold on;plot([400 600],[1 1]*Yfinal(2,2),'--');hold off
axis([0 600 -.02 .02])
title('Step response to 1/6 thrust perturbation')
subplot(233)
plot(Tt(LL2),Yt(LL2,5,1),'LineWidth',2)
ylabel('\gamma');xlabel('time')
hold on;plot([400 600],[1 1]*Yfinal(5,2),'--');hold off
axis([0 600 0 .1])
subplot(234)
plot(Tt2(LL1),Yt2(LL1,1,1),'LineWidth',2)
ylabel('u');xlabel('time')
subplot(235)
plot(Tt2(LL2),Yt2(LL2,2,1)/U0,'LineWidth',2)
ylabel('\alpha (rad)');xlabel('time')
subplot(236)
plot(Tt2(LL1),Yt2(LL1,5,1),'LineWidth',2)
ylabel('\gamma');xlabel('time')
print -depsc ac_fig3a

!ps2pdf ac_fig1.eps
!ps2pdf ac_fig2.eps
!ps2pdf ac_fig3.eps
!ps2pdf ac_fig3a.eps

!cp ac_fig1.eps.pdf ac_fig1.pdf
!cp ac_fig2.eps.pdf ac_fig2.pdf
!cp ac_fig3.eps.pdf ac_fig3.pdf
!cp ac_fig3a.eps.pdf ac_fig3a.pdf