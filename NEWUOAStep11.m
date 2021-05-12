function [] = NEWUOAStep11()
%NEWUOAStep11 对应NEWUOA算法的第十一步

global rho rho_end
if rho==rho_end
    NEWUOAStep13();%To NEWUOAStep13
else
    NEWUOAStep12();%To NEWUOAStep12
end

end