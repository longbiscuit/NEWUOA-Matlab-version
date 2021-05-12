function [] = NEWUOAStep12()
%NEWUOAStep12 对应NEWUOA算法的第十二步

global rho delta Krho
delta=0.5*rho;
rho=0.1*rho;
Krho=0;
NEWUOAStep2();%To Setp2
end