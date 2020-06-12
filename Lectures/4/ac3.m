%
% B747 lateral dynamics
% Reid 188
clear all
%
% metric
   Yv=-1.61e4;Yp=0;Yr=0;
   Lv=-3.062e5;Lp=-1.076e7;Lr=9.925e6;
   Nv=2.131e5;Np=-1.33e6;Nr=-8.934e6;
   
	g=9.81;theta0=0;S=511;cbar=8.324;b=59.64;
   U0=235.9;
   m=2.83176e6/g;cbar=8.324;rho=0.3045;
	Iyy=.449e8;Ixx=.247e8;Izz=.673e8;Ixz=-.212e7;
   
% nondimensional control derivatives
% reid page 207
%
Cyda=0;Cydr=.1146;
Clda=-1.368e-2;Cldr=6.976e-3;
Cnda=-1.973e-4;Cndr=-.1257;

QdS=1/2*rho*U0^2*S;
Yda=QdS*Cyda;Ydr=QdS*Cydr;Lda=QdS*b*Clda;Ldr=QdS*b*Cldr;
Nda=QdS*b*Cnda;Ndr=QdS*b*Cndr;

Ixxp=(Ixx*Izz-Ixz^2)/Izz;
Izzp=(Ixx*Izz-Ixz^2)/Ixx;
Ixzp=Ixz/(Ixx*Izz-Ixz^2);

% reid 289
% full lateral derivatives
% state x=[v p r phi]
%
Alat=[Yv/m Yp/m (Yr/m-U0) g*cos(theta0);
   (Lv/Ixxp + Ixzp*Nv) (Lp/Ixxp + Ixzp*Np) (Lr/Ixxp + Ixzp*Nr) 0;
   (Ixzp*Lv + Nv/Izzp) (Ixzp*Lp + Np/Izzp) (Ixzp*Lr + Nr/Izzp) 0;
   0 1 tan(theta0) 0];

Blat=[1/m 0 0;0 1/Ixxp Ixzp;0 Ixzp 1/Izzp;0 0 0]*[Yda Ydr;Lda Ldr;Nda Ndr];

% Model from rudder to r
G_dr_r=ss(Alat,Blat(:,2),[0 0 1 0],0);
%
% rudder servo 3.33/(s+3.33)
Hr=tf(3.33,[1 3.33]);
%
% expected control gain
Kr=-1.6;

tf(G_dr_r*Hr)

figure(4);clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultlineMarkerSize',18)
pzmap(sign(Kr)*G_dr_r*Hr); 
hh=get(gcf,'children');hhh=get(hh(1),'children')
set(hhh(1),'MarkerSize',18);set(hhh(1),'LineWidth',2)
set(hhh(2),'MarkerSize',18);set(hhh(2),'LineWidth',2)
axis([-3.75 .25 -2 2])
title('Lateral autopilot: r to rudder','FontSize',12)
print -depsc lat_1;jpdf('lat_1')

figure(1);clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultlineMarkerSize',18)
rlocus(G_dr_r*Hr); % note positive sign used here
sgrid([.25 .5 .75],[.25 .5 1])
axis([-2 2 -2 2])
title('Lateral autopilot: r to \delta_r with k>0 ','FontSize',16)
print -depsc lat_2a;jpdf('lat_2a')

figure(1);clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultlineMarkerSize',18)
rlocus(sign(Kr)*G_dr_r*Hr); % note negative sign used here
rr=rlocus(G_dr_r*Hr,Kr); 
hold on; plot(rr+eps*j,'rd','MarkerFace','r');hold off
hh=get(1,'children');hhh=get(hh(1),'children')
set(hhh(1),'MarkerSize',18);set(hhh(1),'LineWidth',2)
set(hhh(2),'MarkerSize',18);set(hhh(2),'LineWidth',2)
sgrid([.25 .5 .75],[.25 .5 1])
axis([-3.75 .25 -2 2])
title('Lateral autopilot: r to \delta_r with k<0 ','FontSize',16)
print -depsc lat_2;jpdf('lat_2')

% clp response to r_comm
%sys=tf(Nrdr*Kpsi,Drdr+Kpsi*Nrdr);
sys_cl=feedback(Kr*G_dr_r*Hr,1,[1],[1],-1)
t=[0:.1:500]';
rc=ones(size(t));
[r1,t]=impulse(sys_cl,t);

% with washout Kwash=[tau 0]/[tau 1]
tau=1/.24; %sec
Nwash=[tau 0];Dwash=[tau 1];
w=logspace(-2,2,250);
Hw=tf(Nwash,Dwash)
Gw=freqresp(Hw,w);
loglog(w,abs(squeeze(Gw)),'LineWidth',2)
axis([.01 20 .01 2])
title('Washout filter with \tau=4.2')
ylabel('|H_w(s)|')
xlabel('Freq (rad/sec)')
print -depsc wash1;jpdf('wash1')

%sysw=tf(conv(Dwash,Nrdr)*Kpsi,conv(Dwash,Drdr)+Kpsi*conv(Nwash,Nrdr));
sys_cl_w=feedback(Kr*G_dr_r*Hr,Hw,[1],[1],-1);
[r2,t]=impulse(sys_cl_w,t);
[u1,t]=lsim(Hw,r2,t);

figure(2);clf;
subplot(211)
plot(t,[r1 r2 ]);axis([0 30 -.5 1 ])
setlines
ylabel('Response r');xlabel('Time')
legend('Without Washout','With Washout')
subplot(212)
plot(t,[Kr*(0-r1) Kr*(0-u1)]);axis([0 30 -.5 1 ])
setlines
ylabel('Control \delta_r');xlabel('Time')
legend('Without Washout','With Washout')
% control signal quickly goes to zero in both cases
% but only one controller has returned r to zero
% the one with a washout produces a zero signal
% even though there is a steady state error in r
print -depsc wash2;jpdf('wash2')

figure(3);clf
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultlineMarkerSize',18)
rlocus(sign(Kr)*G_dr_r*Hr*Hw); % note negative sign used here
rrwash=rlocus(G_dr_r*Hr*Hw,Kr); 
hold on; plot(rr+eps*j,'rd','MarkerFace','r');hold off
hold on; plot(rrwash+eps*j,'m<','MarkerFace','m');hold off
hh=get(1,'children');hhh=get(hh(1),'children')
set(hhh(1),'MarkerSize',18);set(hhh(1),'LineWidth',2)
%set(hhh(2),'MarkerSize',18);set(hhh(2),'LineWidth',2)
sgrid([.25 .5 .75],[.25 .5 1])
axis([-2 1 -1.5 1.5])
title('Lateral autopilot: r to rudder WITH washout filter')
print -depsc wash3;jpdf('wash3')

% roll control
Ha=tf(1,[.15 1]);
G_da_phi=tf(Lda,conv([1 0],[Ixxp -Lp]));
PDr=tf([1/1.5 1],[1]); % locate the PD zero
Kphi=-30;
figure(9);
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultlineMarkerSize',18)
rlocus(-PDr*Ha*G_da_phi)
rrphi=rlocus(PDr*Ha*G_da_phi,Kphi); 
hold on; plot(rrphi+eps*j,'rd','MarkerFace','r');hold off
print -depsc roll1;jpdf('roll1')

%return
%
%
% roll controller reid 289
%Kphi=-30;
%Kr=-1.6;
%
%start with Alat and Blat
%add Psi: dot Psi = r sec(theta)
%state x=[v p r phi Psi]
Alat2=[Alat zeros(4,1);[0 0 sec(theta0) 0 0]];
Blat2=[Blat;[0 0]];
Clat2=eye(5);Dlat2=zeros(5,2);
syslat2=ss(Alat2,Blat2,Clat2,Dlat2);

% from R108, x_e = U0 t 
% from R108, dot y_e = U0 psi cos(theta0) + v
% x=[v p r phi Psi ye]
na=size(Alat2,1);
Alsimul=[Alat2 zeros(na,1);[1 0 0 0 U0*cos(theta0) 0]];
Blsimul=[Blat2;[0 0]];
Clsimul=[Clat2 zeros(5,1);[zeros(1,na) 1]]; % pull out ye as well
Dlsimul=zeros(6,2);

%
% add dynamics on the inputs
%
% aerilon: [0 -1]/[tau_a 1];
tau_a=0.15;
JNa=[0 1];JDa=[tau_a 1];
sysa=tf(JNa,JDa);
% rudder: [0 1]/[tau_r 1];
tau_r=0.3;
JNr=[0 1];JDr=[tau_r 1];
sysr=tf(Kr*JNr,JDr);
sysd=append(sysa,sysr);
sysroll=series(sysd,syslat2);
%
%add the washout filter
%
tau_w=4; %sec
Nwash=[tau_w 0];Dwash=[tau_w 1];
sys_wash=tf(Nwash,Dwash);
sysrollw = feedback(sysroll,sys_wash,2,[3],-1); 
% 

figure(19);
set(gcf,'DefaultLineLineWidth',2)
set(gcf,'DefaultlineMarkerSize',18)
rlocus(-PDr*tf(sysrollw(4,1))) % looks the same as the roll RL above
%rlocus(-PDr*Ha*G_da_phi)
%rrphi=rlocus(PDr*Ha*G_da_phi,Kphi); 
rrphi2=rlocus(PDr*tf(sysrollw(4,1)),Kphi); 
hold on; plot(rrphi2+eps*j,'rd','MarkerFace','r');hold off
print -depsc roll1a;jpdf('roll1a')

syscl = feedback(sysrollw,Kphi*[1/1.5 1],1,[2 4],-1); 
tf(syscl(5,1))
syscl=series(append(tf(Kphi),tf(1)),syscl);% gain on the command input
tf(syscl(5,1))
%return

t=[0:.1:10]';
[yroll,t]=initial(syscl,[0 0 0 15*pi/180 0 0 0 0]',t);
fignum=1;
figure(fignum);clf
set(gcf,'DefaultLineLineWidth',2)
plot(t,yroll(:,[4 2]))
setlines(2)
legend('\phi','p')
title('Reponse to initial roll of 15 degs','FontSize',12)
print -depsc phicomm1;jpdf('phicomm1')
fignum=fignum+1;figure(fignum);clf
set(gcf,'DefaultLineLineWidth',2)
plot(t,yroll(:,[5 3]),t,yroll(:,1)/U0,':')
setlines(2)
legend('\psi','r','\beta')
title('Reponse to initial roll of 15 degs','FontSize',12)
print -depsc phicomm1a;jpdf('phicomm1a')

t=[0:.1:30]';
fact=15*pi/180;
[ystep,t]=step(syscl,t);
fignum=fignum+1;figure(fignum);clf
set(gcf,'DefaultLineLineWidth',2)
plot(t,ystep(:,[4 2])*fact)
xlabel('Time (sec)')
setlines(2)
legend('\phi','p')
title('Reponse to step input \phi_c=15 degs (0.26 rad)','FontSize',12)
print -depsc phicomm2;jpdf('phicomm2')
fignum=fignum+1;figure(fignum);clf
set(gcf,'DefaultLineLineWidth',2)
plot(t,ystep(:,[5])*fact,t,ystep(:,[3])*2,t,ystep(:,1)/U0*fact,':')
xlabel('Time (sec)')
setlines(2)
legend('\psi','r','\beta')
title('Reponse to step input \phi_c=15 degs (0.26 rad)','FontSize',12)
print -depsc phicomm2a;jpdf('phicomm2a')

%add the heading controller
%phi_c = (U0/g/tau_h) * (Psi_d - Psi);
% input is now going to be Psi_d
%
tau_h=12; % time constant to follow the heading command
rgain=(U0/g/tau_h);
fignum=fignum+1;figure(fignum);clf
rlocus(syscl(5,1))
rr=rlocus(syscl(5,1),rgain)
hold on; plot(rr+eps*j,'rd','MarkerFace','r');hold off
axis([-2 2 -1.25 1.25]);title('Heading Command Loop \tau_1=12 sec','FontSize',14)
print -depsc phicomm3;jpdf('phicomm3')

syscl2=feedback(syscl,rgain,[1],[5],-1)
syscl2=series(append(tf(rgain),tf(1)),syscl2);% gain on the command input
t=[0:.1:50]';
fact=15*pi/180; %step of 15 deg in heading
[yhead,t]=step(syscl2(:,1),t);
fignum=fignum+1;figure(fignum);clf
plot(t,yhead(:,[4 2])*fact)
xlabel('Time (sec)')
setlines(2)
legend('\phi','p')
print -depsc phicomm3a;jpdf('phicomm3a')
fignum=fignum+1;figure(fignum);clf
plot(t,yhead(:,[5])*fact,t,yhead(:,[3]),t,yhead(:,1)/U0*fact,'m:',t,fact*ones(size(t)))
xlabel('Time (sec)')
setlines(2)
legend('\psi','r','\beta','\psi_c')
print -depsc phicomm3b;jpdf('phicomm3b')
%

% now determine psi_d from the track that we would like to follow
% 
tau_2 = 30; %sec. this is the timescale for the path follower
%

return

%run ac3_lat
%ac3out

T=ac3out(:,1);
V=ac3out(:,2);
P=ac3out(:,3);
R=ac3out(:,4);
PHI=ac3out(:,5);
PSI=ac3out(:,6);
YREF=ac3out(:,7);
YE=ac3out(:,8);
PSIcomm=ac3out(:,9);
PHIcomm=ac3out(:,10);

figure(1);clf
plot(T,[YREF YE]);
f=get(gcf,'children')
set(f,'Ydir','reverse')
setlines(2)
ylabel('Path Y_e')
legend('Y_{ref}','Y_{a/c}')
print -depsc ac3_lat_path;jpdf('ac3_lat_path')

figure(2);clf
plot(T,[PHI PHIcomm]);
setlines(2)
ylabel('\Phi and \Phi_c')
legend('\Phi','\Phi_c')
print -depsc ac3_lat_phi;jpdf('ac3_lat_phi')

figure(3);clf
plot(T,[PSI PSIcomm]);
setlines(2)
ylabel('\Psi and \Psi_c')
legend('\Psi','\Psi_c')
print -depsc ac3_lat_psi;jpdf('ac3_lat_psi')
