function [] = NEWUOAStep6()
%NEWUOAStep6 对应NEWUOA算法的第六步

global RATIO
if RATIO>=0.1
    NEWUOAStep2();%To NEWUOAStep2
else
    NEWUOAStep7();%To NEWUOAStep7
end

end