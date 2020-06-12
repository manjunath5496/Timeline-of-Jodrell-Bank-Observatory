% newr.m 
% Analyze tracking algorithm by Park et al
% AIAA GNC 2004 
% 
% Assumes that ac3.m has been run to generate syscl
% Jonathan How
% MIT 16.333 Fall 2004
%
%
close all
dt=1; % time step for the simulation
U0=235.9;
path=[];

jcase=1;
% 2 path cases considered
if jcase==1
    t=[0:5*dt:2500]';
    omp=.0025;
    path=24000*[sin(omp*t) (1-cos(omp*t))];
    xe=0;ye=1500;
    X=[.1 0 0 0*pi/180 0*pi/180 0 0 0]';
else
    t=[0:dt:1350]';
    path(:,1)=U0*t;
    omp=.005;
    path(:,2)=500*(-cos(2*pi*omp*t)+1).*exp(.002*t);
    xe=0;ye=-1000;
    X=[.1 0 0 -15*pi/180 -15*pi/180 0 0 0]';
end

% Discretize the dynamics with time step dt
% system has the inner yaw and roll loops closed
[A,B,C,D]=ssdata(syscl);
syscld=c2d(ss(A,B,C,D),dt);
[Ad,Bd,Cd,Dd]=ssdata(syscld);
Bd=Bd(:,1);Dd=Dd(:,1); % only need first input

% bank angle limit
philim=30;
%
%inputs are phi_d and 0
%state x=[v p r phi Psi xx xx xx]
L1=2000; % look ahead distance
store=[];

% find the point on the path L1 ahead
ii=find((xe-path(:,1)).^2+(ye-path(:,2)).^2 < L1^2);
iii=max(ii);
%
%
%
kk=1;
while (~isempty(iii)) & (kk< length(t)-1)
    kk=kk+1;    
    aim_point=path(iii,:);
    
    xedot=U0*cos(X(5))-X(1)*cos(X(4))*sin(X(5));
    yedot=U0*sin(X(5))+X(1)*cos(X(4))*cos(X(5));
    
    v1=[xedot yedot]';
    v2=[aim_point(1)-xe aim_point(2)-ye]';
    v1=v1/norm(v1);
    v2=v2/norm(v2);
    [v1 v2];
    temp=cross([v1;0],[v2;0]);
    eta=acos(v1'*v2)*sign(temp(3));
    phi_d=atan(2*U0^2/L1*sin(eta)/9.81);
    phi_d=max(min(phi_d,philim*pi/180),-philim*pi/180);
    
    store=[store;[t(kk) X' xe ye phi_d v2']];
    % propagate forward a step
    X=Ad*X+Bd*phi_d;
    xe=xe+xedot*dt;
    ye=ye+yedot*dt;
    
    ii=find((xe-path(:,1)).^2+(ye-path(:,2)).^2 < L1^2);
    iii=max(ii);
end

figure(1);clf
plot(store(:,11),store(:,10),'m');
hold on;plot(path(:,2),path(:,1),'g');
legend('veh','Path');xlabel('Y_e');
ylabel('x_e');setlines(2);hold off
if jcase==1
    axis('square');axis('equal')
    print -depsc park_1; jpdf('park_1')
else
    orient tall
    print -depsc park_1a; jpdf('park_1a')
end

figure(2);clf
plot(store(:,1),store(:,[5 12])*180/pi);
axis([0 t(kk) -philim*1.1 philim*1.1])
xlabel('time');ylabel('\phi and \phi_d');
legend('\phi','\phi_d');setlines(2)
if jcase==1
    print -depsc park_2; jpdf('park_2')
else
    print -depsc park_2a; jpdf('park_2a')
end
return
