function [] = NEWUOAStep15()
%NEWUOAStep15 对应NEWUOA算法的第十五步

global RATIO delta rho
delta=0.1*delta;
RATIO=-1;
if delta<=1.5*rho
    delta=rho;
end
NEWUOAStep7();%To Setp7
end