function [] = NEWUOAStep4()
%NEWUOAStep4 对应NEWUOA算法的第四步
% 计算比例系数RATIO，以此作为迭代delta的依据。
% 把MOVE设置为需要被替换的点
global m n RATIO Fn Xn F F_times delta opt d NORMD rho MOVE Qnew H x0 w Hw beta Fnew Krho
Fopt=Fn(opt);
xopt=Xn(:,opt);
xnew=xopt+d;%新生成的点
Fnew=F(xnew);
F_times=F_times+1;
Krho=Krho+1;
%%
%计算RATIO
dQ=Fopt-Qnew;
dF=Fopt-Fnew;
if dQ<=0 %%二次模型没有下降
    if dF>0 %但是实际值是下降的
        RATIO=0;%去修改模型
    else%实际值也没有下降
        RATIO=-1;
    end
else
    RATIO=dF/dQ;
end
%%
%更新delta
if RATIO<=0.1
    deltaint=0.5*NORMD;
else
    if RATIO<=0.7
        deltaint=max([NORMD,0.5*delta]);
    else
        deltaint=max([2*NORMD,0.5*delta]);
    end
end
if deltaint>1.5*rho
    delta=deltaint;
else
    delta=rho;
end
%%
%挑选出被替代的点
    T=1:m;%待挑选的点 
if dF>0
    Case=1;
    xstar=xopt+d;%之后的最优值点
else
    xstar=xopt;%之后的最优值点
    T(opt)=[];
    Case=-1;
end
w=zeros(m+n+1,1);
dxx0=xnew-x0;
for i=1:m
    w(i)=0.5*(((Xn(:,i)-x0)'*(dxx0))^2);
end
w(m+1)=1;
w(m+2:m+n+1)=dxx0;
Hw=H*w;
beta=0.5*((dxx0'*dxx0)^2)-w'*Hw;
Sigma=zeros(m,1);
Weight=zeros(m,1);
M1=max([0.1*delta,rho]);
if Case>0
    for i=1:m
        alpha=H(i,i);
        tau=Hw(i);
        Sigma(i)=alpha*beta+tau*tau;
        Weight(i)=max([1,(norm(Xn(:,i)-xstar)/M1)^6]);
    end
    [~,tstar]=max(Weight.*abs(Sigma));
    MOVE=tstar;
else
    for i=1:m-1
        alpha=H(T(i),T(i));
        tau=Hw(T(i));
        Sigma(T(i))=alpha*beta+tau*tau;
        Weight(T(i))=max([1,(norm(Xn(:,T(i))-xstar)/M1)^6]);
    end
    [M2,tstar]=max(Weight.*abs(Sigma));
    if M2<=1
        MOVE=0;
    else
        MOVE=tstar;
    end
end



%%
 NEWUOAStep5();%To Setp5







end