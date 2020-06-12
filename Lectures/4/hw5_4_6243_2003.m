function hw5_4_6243_2003(n)
% function hw5_4_6243_2003(n)
%
% calculates numerically the evolution matrix M(2*pi)
% of system dx/dt=A(t)x(t), where A(t)=[0 1;-3sin^2(t) -1]
% using n step integration

M=eye(2);
T=pi;
for k=1:n,
    M=expm([0 1;-3*sin(k*T/n) -1]*(T/n))*M;
end
M
eig(M)
det(M)
exp(-pi)