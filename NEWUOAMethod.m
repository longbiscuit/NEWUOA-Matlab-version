function [Fopt,xopt] = NEWUOAMethod(Fobj,M,N,xbeg,rhobeg,rhoend,Max)
%NEWUOAMETHOD NEWUOA算法的函数
%   F表示目标函数,m是插值点个数,n是问题维数,rho_beg是初始的rho,rho_end是终止的rho,Max是迭代的最大次数
global Xn Fn n m F_times rho_beg rho_end x0 opt xb c g Gamma gamma H F rho delta Krho D3 QF3 CRVMIN d NORMD Qnew 
global RATIO MOVE w Hw beta Fnew DIST XXopt NXX 
F=Fobj;
m=M;
n=N;
xb=xbeg;
rho_beg=rhobeg;
rho_end=rhoend;

NEWUOAStep1();
Fopt=Fn(opt);
xopt=Xn(:,opt);
end

