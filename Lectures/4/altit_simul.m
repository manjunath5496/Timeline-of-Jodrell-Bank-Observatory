%
% altitude hold controller
%
clear all
prt=1;
for ii=1:9
  figure(ii);clf;
  set(gcf,'DefaultLineLineWidth',2);
  set(gcf,'DefaultlineMarkerSize',10)
end

Xu=-1.982e3;Xw=4.025e3;
Zu=-2.595e4;Zw=-9.030e4;Zq=-4.524e5;Zwd=1.909e3;
Mu=1.593e4;Mw=-1.563e5;Mq=-1.521e7;Mwd=-1.702e4;

g=9.81;theta0=0;S=511;cbar=8.324;
U0=235.9;Iyy=.449e8;m=2.83176e6/g;cbar=8.324;rho=0.3045;
Xdp=.3*m*g;Zdp=0;Mdp=0;
Xde=-3.818e-6*(1/2*rho*U0^2*S);Zde=-0.3648*(1/2*rho*U0^2*S);
Mde=-1.444*(1/2*rho*U0^2*S*cbar);;
%
%
% use the full longitudinal model 
% 
% x=[u w q theta h];
% u=[de;dt];
sen_u=1;sen_w=2;sen_q=3;sen_t=4;sen_h=5;sen_de=6;sen_dt=7;
act_e=1;act_t=2;
%
% dot h = U0 (theta - alpha) = U0 theta - w
%
A=[Xu/m Xw/m 0 -g*cos(theta0) 0;[Zu Zw Zq+m*U0 -m*g*sin(theta0)]/(m-Zwd) 0;
[Mu+Zu*Mwd/(m-Zwd) Mw+Zw*Mwd/(m-Zwd) Mq+(Zq+m*U0)*Mwd/(m-Zwd) ...
     -m*g*sin(theta0)*Mwd/(m-Zwd) 0]/Iyy;
[ 0 0 1 0 0];[0 -1 0 U0 0]]; 
B=[Xde/m Xdp/m;Zde/(m-Zwd) Zdp/(m-Zwd);(Mde+Zde*Mwd/(m-Zwd))/Iyy ...
  (Mdp+Zdp*Mwd/(m-Zwd))/Iyy;0 0;0 0];
C=[eye(5);zeros(2,5)];
D=[zeros(5,2);[eye(2)]]; %last 2 outputs are the controls
%
% add actuator dynamics
%
% elevator
%tau_e=.1;%sec
tau_e=.25;%sec
JNe=1;JDe=[tau_e 1];
syse=tf(JNe,JDe);
%
% thrust 
%
tau_t=3.5;%sec
JNt=1;JDt=[tau_t 1];
syst=tf(JNt,JDt);
sysd=append(syse,syst);
syslong=ss(A,B,C,D);
syslong2=series(sysd,syslong);
[A,B,C,D]=ssdata(syslong2);
na=size(A);
%
% q and theta loops
%
% inner loop on theta/q
% u=kq*q+Kth*theta + de_c
K_th=1;K_q=1.95*K_th;
if 1
   figure(1);clf;%orient tall
   sgrid([.5 .707]',[1]);
   rlocus(tf(-[1.95 1],1)*tf(syslong2(sen_t,act_e)))
   r_th=rlocus(tf(-[1.95 1],1)*tf(syslong2(sen_t,act_e)),K_th)
   hold on
   plot(r_th+eps*j,'bd','MarkerFace','b');axis([-3,1,-3,3]);hold off
   axis([-3 .1 -3 3])
   title('with q and theta FB to \delta_e');
   if prt
      print -depsc alt_sim1.eps
      jpdf('alt_sim1')
   end
end
% close inner loop
syscl=feedback(syslong2,[K_q K_th],act_e,[sen_q sen_t],1)
[eig(syscl) [nan;r_th;nan]]
[Acl,Bcl,Ccl,Dcl]=ssdata(syscl);
%
% engine loops
% mostly the phugoid ==> design on top of the short period
%
% system with q/th loop feedback
figure(9);clf
% interested in vel output and engine input
sys_u=syscl(sen_u,act_t);
Kun=[1/0.2857 1]; % comp zero on plant pole
Kud=[1/(0.2857*5) 1];
[Au,Bu,Cu,Du]=tf2ss(Kun,Kud);
%K_u=.075;
%K_u=.35;
K_u=.1;
rlocus(sys_u*tf(Kun,Kud));
rr_u=rlocus(sys_u*tf(Kun,Kud),K_u);
hold on;plot(rr_u+eps*j,'bd','MarkerFace','b');hold off;
 
if 1
   na=size(Acl,1);
   Aclt=[Acl Bcl(:,act_t)*Cu;zeros(1,na)  Au];
   Bclt=[[Bcl(:,act_e);0] [Bcl(:,act_t)*Du;Bu]];
   Cclt=[Ccl zeros(size(Ccl,1),1)];
   Dclt=[zeros(size(Ccl,1),2)];
   figure(2);clf;axis([-.5 .05,-.4,.4]);sgrid([.5 .707]',[.05]);hold on;
   rlocus(Aclt,Bclt(:,act_t),Cclt(sen_u,:),Dclt(sen_u,act_t));
   axis([-2 .05,-2,2])
   r_u=rlocus(Aclt,Bclt(:,act_t),Cclt(sen_u,:),Dclt(sen_u,act_t),K_u)';
   [[r_th;NaN;NaN;NaN] r_u]
   hold on;plot(r_u+eps*j,'bd','MarkerFace','b'),hold off
   title('with u FB to \delta_t')
   if prt
      print -depsc alt_sim2.eps
     jpdf('alt_sim2')
   end
end
Acl2=Aclt-Bclt(:,act_t)*K_u*Cclt(sen_u,:);
Bcl2=Bclt;
Ccl2=Cclt;
Dcl2=Dclt;
sysclt=ss(Acl2,Bcl2,Ccl2,Dcl2);
[eig(sysclt) rr_u]
%
% now close loop on altitude
%
% de_c = kh*(h_c-h)
if 1
    % design lead by canceling troubling plant pole 
    % zero located p*8
   tt=eig(sysclt);[ee,ii]=min(abs(tt+.165));
   Khn=[1/abs(tt(ii)) 1];Khd=[1/(8*abs(tt(ii))) 1];
   K_h=-1*.00116;
   Gc_eng=tf(Khn,Khd);
   Loopt=series(append(Gc_eng,tf(1,1)),sysclt);    
   figure(3);clf;
   sgrid([.5 .707]',[.1:.1:1]);
   hold on;
   rlocus(sign(K_h)*Loopt(sen_h,act_e));
   axis([-1 .1, -2 2])
   r_h=rlocus(Loopt(sen_h,act_e),K_h)
   hold on;plot(r_h+eps*j,'bd','MarkerFace','b'),hold off
   title('with h FB to \delta_e command')
   if prt
      print -depsc alt_sim22.eps
     jpdf('alt_sim22')
   end
end
syscl3=feedback(series(append(tf(K_h,1),tf(1,1)),Loopt),[1],act_e,sen_h,-1);
[r_h eig(syscl3)]
%
% try using the thrust act to control altitude
%
if 1
%   K_hu=.154;
   K_hu=0.025;
   tt=eig(sysclt);[ee,ii]=min(abs(tt+.165));
   Khnu=[1/abs(tt(ii)) 1];Khdu=[1/(5*abs(tt(ii))) 1];
   Gc_engu=tf(Khnu,Khdu);
   Loopt=series(append(tf(1,1),Gc_engu),sysclt);    
   figure(7);clf;axis([-.6 .1,-.5,.5]);sgrid([.5 .707]',[.05]);hold on;
   rlocus(sign(K_hu)*Loopt(sen_h,act_t));
   axis([-1 .1,-2,2])
   r_hu=rlocus(Loopt(sen_h,act_t),K_hu)
   hold on;plot(r_hu+eps*j,'bd','MarkerFace','b'),hold off
   title('with h FB to \delta_t command')
   if prt
      print -depsc alt_sim22a.eps
     jpdf('alt_sim22a')
   end
end
syscl4=feedback(series(append(tf(1,1),tf(K_hu,1)),Loopt),[1],act_t,sen_h,-1);
[r_hu eig(syscl4)]

t=[0:.01:18]';
[ystep,t]=initial(syscl3,[0 0 0 0 90 0 0 0 0],t);
figure(4);clf;orient tall;
subplot(211)
ll=1:length(t);
U=ystep(:,sen_u);W=ystep(:,sen_w);q=ystep(:,sen_q);
TH=ystep(:,sen_t);H=ystep(:,sen_h);
DELe=ystep(:,sen_de);DELt=ystep(:,sen_dt);
plot(t(ll),H(ll),t(ll),5*U(ll),t(ll),TH(ll)*180/pi);setlines(2)
title('Initial Condition response to 90m altitude error. E-FB')
legend('H','5*U','\theta deg');grid
subplot(212)
plot(t(ll),DELe(ll)*180/pi,t(ll),20*DELt(ll));
setlines(2)
legend('\delta_e input','20*\delta_t input');grid
axis([0 18 -20 20])
ylabel('Commands (degs)')
title('Initial Condition response to 90m altitude error')
if prt
   print -depsc alt_sim3.eps
   jpdf('alt_sim3')
end
%
t=[0:.01:18]';
[ystep,t]=initial(syscl4,[0 0 0 0 90 0 0 0 0],t);
figure(14);clf;orient tall;
subplot(211)
ll=1:length(t);
U=ystep(:,sen_u);W=ystep(:,sen_w);q=ystep(:,sen_q);
TH=ystep(:,sen_t);H=ystep(:,sen_h);
DELe=ystep(:,sen_de);DELt=ystep(:,sen_dt);
plot(t(ll),H(ll),t(ll),U(ll),t(ll),TH(ll)*180/pi);setlines(2)
title('Initial Condition response to 90m altitude error. U - FB')
legend('H','U','\theta deg');grid
subplot(212)
plot(t(ll),DELe(ll)*180/pi,t(ll),DELt(ll));
setlines(2)
legend('\delta_e input','\delta_t input');grid
axis([0 18 -20 20])
ylabel('Command inputs (degs)')
if prt
   print -depsc alt_sim3a.eps
   jpdf('alt_sim3a')
end

figure(5);clf
%Path to follow
Hmax=1500;
Hc=[zeros(10,1);Hmax/100*[0:1:100]';Hmax*ones(50,1);Hmax/100*[100:-1:0]';zeros(50,1)];
T=0:length(Hc)-1;
[Yh,T]=lsim(syscl3(:,1),Hc,T);
U=Yh(:,sen_u);W=Yh(:,sen_w);q=Yh(:,sen_q);TH=Yh(:,sen_t);H=Yh(:,sen_h);
DELe=Yh(:,sen_de);DELt=Yh(:,sen_dt);
plot(T,[Hc],T,H,'--');setlines(2)
axis([0 300 -100 1600])
legend('h_c','h')
title('Altitude controller: elevator FB')
ylabel('height')
xlabel('time')
if prt
   print -depsc alt_sim4.eps
   jpdf('alt_sim4')
end

figure(6);clf;
subplot(211)
ll=1:length(T);
plot(T(ll),U(ll),T(ll),W(ll)/U0*180/pi,T(ll),TH(ll)*180/pi);setlines(2)
title('Linear Response - elevator-FB')
legend('U','\alpha deg','\theta deg');grid
axis([0 300 -5 5])
subplot(212)
plot(T(ll),DELe(ll)*180/pi,T(ll),10*DELt(ll));
setlines(2)
axis([0 300 -3 3])
legend('\delta_e','10*\delta_t');grid
%axis([0 250 -4 4])
if prt
   print -depsc alt_sim5.eps
   jpdf('alt_sim5')
end

figure(15);clf
%Path to follow
Hmax=1500;
Hc=[zeros(10,1);Hmax/100*[0:1:100]';Hmax*ones(50,1);Hmax/100*[100:-1:0]';zeros(50,1)];
T=0:length(Hc)-1;
[Yh,T]=lsim(syscl4(:,2),Hc,T);
U=Yh(:,sen_u);W=Yh(:,sen_w);q=Yh(:,sen_q);TH=Yh(:,sen_t);H=Yh(:,sen_h);
DELe=Yh(:,sen_de);DELt=Yh(:,sen_dt);
plot(T,[Hc],T,H,'--');setlines(2)
axis([0 300 -100 1600])
legend('h_c','h')
title('Altitude controller: thrust FB')
ylabel('height')
xlabel('time')
if prt
   print -depsc alt_sim4a.eps
   jpdf('alt_sim4a')
end

figure(16);clf;
subplot(211)
ll=1:length(T);
plot(T(ll),U(ll),T(ll),5*W(ll)/U0*180/pi,T(ll),5*TH(ll)*180/pi);setlines(2)
title('Linear Response - thrust-FB')
legend('U','5*\alpha deg','5*\theta deg');grid
axis([0 300 -40 40])
subplot(212)
plot(T(ll),DELe(ll)*180/pi,T(ll),DELt(ll));
setlines(2)
legend('\delta_e','\delta_t');grid
axis([0 300 -2 2])
%axis([0 250 -4 4])
if prt
   print -depsc alt_sim5a.eps
   jpdf('alt_sim5a')
end

return

sT=alt(:,1);u=alt(:,2);w=alt(:,3);
q=alt(:,4);theta=alt(:,5);h=alt(:,6);
dele=alt(:,7);delt=alt(:,8);href=alt(:,9);

figure(17);clf;orient tall
ll=1:length(sT)-2;
subplot(311)
title('Simulink hard path following: elevator FB')
plot(sT(ll),u(ll),sT(ll),w(ll)/U0*180/pi,sT(ll),theta(ll)*180/pi);setlines(2)
legend('U','\alpha deg','\theta deg');grid
subplot(312)
plot(sT(ll),[href(ll) h(ll)])
setlines(2)
legend('H_{ref}','H');grid
axis([0 300 -100 1600])
subplot(313)
plot(sT(ll),dele(ll)*180/pi,sT(ll),10*delt(ll));
setlines(2)
legend('\delta_e','10*\delta_t');grid
%axis([0 60 -1000 600])
if prt
   print -depsc alt_sim6.eps
   jpdf('alt_sim6')
end

sTu=altu(:,1);uu=altu(:,2);wu=altu(:,3);
qu=altu(:,4);thetau=altu(:,5);hu=altu(:,6);
deleu=altu(:,7);deltu=altu(:,8);hrefu=altu(:,9);

figure(18);clf;orient tall
ll=1:length(sTu)-2;
subplot(311)
title('Simulink hard path following: thruster FB')
plot(sTu(ll),uu(ll),sTu(ll),wu(ll)/U0*180/pi,sTu(ll),thetau(ll)*180/pi);setlines(2)
legend('U','\alpha deg','\theta deg');grid
subplot(312)
plot(sTu(ll),[hrefu(ll) hu(ll)])
setlines(2)
legend('H_{ref}','H');grid
axis([0 300 -100 1600])
subplot(313)
plot(sTu(ll),deleu(ll)*180/pi,sTu(ll),deltu(ll));
setlines(2)
legend('\delta_e','\delta_t');grid
%axis([0 60 -1000 600])
if prt
   print -depsc alt_sim6a.eps
   jpdf('alt_sim6a')
end
