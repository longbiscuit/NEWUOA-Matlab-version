function [] = NEWUOAStep10()
%NEWUOAStep10 对应NEWUOA算法的第十步

global delta RATIO NORMD rho
if max(NORMD,delta)<=rho&&RATIO>=0
    NEWUOAStep11();%To NEWUOAStep11
else
    NEWUOAStep2();%To NEWUOAStep2
end

end