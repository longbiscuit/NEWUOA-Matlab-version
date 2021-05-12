function [] = NEWUOAStep3()
%NEWUOAStep3 对应NEWUOA算法的第三步

global NORMD rho
if NORMD>0.5*rho
    NEWUOAStep4();%To NEWUOAStep4
else
    NEWUOAStep14();%To NEWUOAStep14
end

end

