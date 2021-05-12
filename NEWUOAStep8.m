function [] = NEWUOAStep8()
%NEWUOAStep8 对应NEWUOA算法的第八步

global DIST delta
if DIST>=2*delta
    NEWUOAStep9();%To NEWUOAStep9
else
    NEWUOAStep10();%To NEWUOAStep10
end

end