%Calculate the Static Stability Params for the 
%Navion General Aviation A/C
% 16.333
% Adapted from Nelson, page 57
%clear all

W=2750 %lb
V=176 % ft/sec
X_cg=0.295 %cbar
h=0.295
S=184 % ft^2
b=33.4 % ft
cbar = 5.7 % ft
S_t=43 %ft^2
l_t=16 %ft
h_t=l_t/cbar+h

% Wing
C_m_ac_w=-0.116
C_l_alp_w=0.097  %/deg
alp_0_L=-5 % deg
X_ac=0.25 %cbar
h_nbar=0.25
i_w=1 %deg

%Tail
C_l_alp_t=0.088 %/deg
C_m_ac_t=0
%i_t=-1% deg
eta=1;
b_t=2.2/5.5*b;

V_H=l_t*S_t/cbar/S
AR=b^2/S
AR_t=b_t^2/S_t
d2r=180/pi;

% convert wing section lifts to 3D wing
C_L_alp_w=C_l_alp_w*d2r/(1+C_l_alp_w*d2r/(pi*AR))
C_L_alp_t=C_l_alp_t*d2r/(1+C_l_alp_t*d2r/(pi*AR_t))

% Wing contribution
% from 2-4
% C_m_cg_w=C_L_w(h-h_nbar)+C_m_ac_w;
%
% with C_L_w=C_L_0_w + C_L_alp_w alp_w
% where C_L_0_w = lift coeff at zero angle of attack
C_L_0_w = C_L_alp_w * abs(alp_0_L/d2r)
%
C_m_cg_0_w=C_L_0_w*(h-h_nbar)+C_m_ac_w
C_m_cg_alp_w=C_L_alp_w*(h-h_nbar)

% for downwash
% eps = eps_0 + eps_alp alp_w
eps_0=2*C_L_0_w/pi/AR
eps_alp=2*C_L_alp_w/pi/AR

% Tail from 2-8
C_m_cg_0_t=eta*V_H*C_L_alp_t*(eps_0+i_w/d2r-i_t/d2r)
C_m_cg_alp_t=-eta*V_H*C_L_alp_t*(1-eps_alp)

% total
C_m_cg_0=C_m_cg_0_w+C_m_cg_0_t
C_m_cg_alp=C_m_cg_alp_w+C_m_cg_alp_t

figure(1);clf
alp=[0:1:12]';
plot(alp,C_m_cg_0_w+C_m_cg_alp_w*alp/d2r,'LineWidth',2)
hold on
plot(alp,C_m_cg_0_t+C_m_cg_alp_t*alp/d2r,'m--','LineWidth',2)
plot(alp,C_m_cg_0+C_m_cg_alp*alp/d2r,'r:','LineWidth',2)
plot([0 12],[0 0],'k','LineWidth',2)
legend('Wing','Tail','Total')
%hold off
xlabel('\alpha (deg)')
ylabel('C_{m_{cg}}')

% neutral point
C_L_alp_T=C_L_alp_w+eta*S_t/S*C_L_alp_t*(1-eps_alp)
gamma=C_L_alp_T/C_L_alp_w
hnp=(h_nbar+(gamma-1)*h_t)/gamma

return

clear all
Hnp=[];kk=0;CMT=[];CMt=[];CMw=[];
for i_t=[-2:1]
    kk=kk+1;
    navion_1
    Hnp=[Hnp hnp]
    CMw=[CMw C_m_cg_0_w+C_m_cg_alp_w*alp/d2r];
    CMt=[CMt C_m_cg_0_t+C_m_cg_alp_t*alp/d2r];
    CMT=[CMT C_m_cg_0+C_m_cg_alp*alp/d2r];
    hold on
end

figure(2);clf
plot(alp,CMw,'LineWidth',2)
hold on
plot(alp,CMt,'--','LineWidth',2)
plot(alp,CMT,'r:','LineWidth',2)
legend('Wing','Tail','Total')
plot([0 12],[0 0],'k','LineWidth',2)
hold off
xlabel('\alpha (deg)')
ylabel('C_{m_{cg}}')
