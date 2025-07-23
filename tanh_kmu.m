function [material]=tanh_kmu(z_end,kmax,kmin,rv,N)

% from deng et al. 2022
% make the distribution of K, MU as Eq.9 of Goriely and Tabor 2003


L = 1;
A=kmax-kmin;
beta=kmin;

% RX - THESE PROPERTIES MIGHT CHANGE BASED ON Z_END, THESE ARE FOR Z_END=4

sigma1=L*(z_end-1);
plane=0.4;%how wide the constant part we want
alpha=0.5*(1-2*plane)*2;
% material property distribution
M=@(z)0.5*A*(1-tanh((z-sigma1)/alpha))+beta;

material=ones(N,1);
for i=1:N
    material(i)=M(0.5*(rv(1,i)+rv(1,i+1)));
end

