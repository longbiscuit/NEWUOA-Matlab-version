clear all
%%
%测试函数
F_test=@(x)(x^2+100*((x^2-1)^2));
%%
%初始化参数
n=1;
m=2*n+1;
rho_beg=1;
rho_end=10^(-6);
Max=10;
xb=100;
%%
global Xn Fn n m F_times rho_beg rho_end x0 opt xb c g Gamma gamma H F rho delta Krho D3 QF3 CRVMIN d NORMD Qnew 
global RATIO MOVE w Hw beta Fnew DIST XXopt NXX 
[Fopt,xopt] = NEWUOAMethod(F_test,m,n,xb,rho_beg,rho_end,Max);