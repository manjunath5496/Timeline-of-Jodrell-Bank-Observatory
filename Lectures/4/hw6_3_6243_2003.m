function hw6_3_6243_2003(x0,x1)
% function hw6_3_6243_2003(x0,x1)
%
% moving the state of system
%    x1'=cos(x3)u1
%    x2'=sin(x3)u1
%    x3'=x4*u1
%    x4'=u2
% from x0 to x1 (default: randomly generated)

if nargin<2, x1=0.5*randn(4,1); end
if nargin<1, x0=0.5*randn(4,1); end

dt=0.01;   % simulation step
tol=0.01;    % positioning accuracy
N=4000;    % maximal number of steps
x=zeros(4,N);  % to keep the state values
u=zeros(2,N);  % to keep the control values
kk=zeros(1,N); % to keep strategies
t=1;        % current time
x(:,t)=x0(:);
x1=x1(:);

while t<N,
    % g=[g1 g2 g3 g4]
    g=[cos(x(3,t)) 0 0 -sin(x(3,t));sin(x(3,t)) 0 0 cos(x(3,t));x(4,t) 0 1 0;0 1 0 0];
    v=g'*(x1-x(:,t))./[(1+x(4,t)^2);1;1;1];  % normalzation
    if any(abs(v)>tol),
        [y,k]=max(abs(v));
        du=S(k,v(k));
    else
        du=zeros(2,N-t);
    end
    for t1=1:min(size(du,2),N-t),     % simulation
        u(:,t+1)=du(:,t1);
        x(:,t+1)=x(:,t)+dt*[cos(x(3,t)) 0;sin(x(3,t)) 0;x(4,t) 0;0 1]*u(:,t+1);
        kk(t+1)=k;
        t=t+1;
    end
end

tt=dt*(0:N-1);
close(gcf)
for k=1:4,
    subplot(5,1,k);
    plot(tt,x(k,:),tt,x1(k)); grid
end
subplot(5,1,5);plot(tt,0.9*kk);grid;


function u=S(k,v)
% control needed to implement the differential flow for m steps along gk
dt=0.01;
switch k,
    case 1,
        m=sign(v)*(floor(abs(v)/dt)+1);
        u=[repmat(sign(m),1,abs(m));zeros(1,abs(m))];
    case 2,
        m=sign(v)*(floor(abs(v)/dt)+1);
        u=[zeros(1,abs(m));repmat(sign(m),1,abs(m))];        
    case 3,    
        v1=sign(v)*sqrt(abs(v));
        u=[S(2,v1) S(1,abs(v1)) S(2,-v1) S(1,-abs(v1))];
    case 4,
        v1=sign(v)*sqrt(abs(v));
        u=[S(3,v1) S(1,abs(v1)) S(3,-v1) S(1,-abs(v1))]; 
end
        
        
            