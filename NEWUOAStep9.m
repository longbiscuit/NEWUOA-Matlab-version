function [] = NEWUOAStep9()
%NEWUOAStep9 对应NEWUOA算法的第九步
%BIGLAGandBIGDENOfNEWUOA 移动插值点的方法
%   目标函数为|l_t(xopt+d)|
%   约束条件为 ||d||<=deltah
%%
global d NORMD w Hw rho delta DIST W H Xn MOVE opt x0 XXopt NXX n m RATIO Fnew F_times Krho
%BIGLAG
deltah=max([min([0.1*DIST,0.5*delta]),rho]);
N=n;
t=MOVE;
xt=Xn(:,t);
xopt=Xn(:,opt);
X=W(m+2:m+n+1,1:m);
%%
%l_t(x)=c+g*(x-x0)+0.5*(x-x0)'*G*(x-x0)
Lambda=diag(H(1:m,t));
c=H(m+1,t);
g=H(m+2:m+n+1,t);
G=X*Lambda*X';
%%
Dl=g+G*(xopt-x0);%Dlt(xopt);
eta=@(u)(H(1:m,t).*(X'*u));
DDlu=@(u)( X*eta(u) );
xx=xt-xopt;
d1=deltah*(xx/norm(xx));
d2=-d1;
DDld=DDlu(d1);
dDDld=d1'*DDld;
dDl=d1'*Dl;
if dDl*dDDld>0 %得到d0
    d=d1;
    lxd=abs(dDl+0.5*dDDld);
else
    d=d2;
    lxd=abs(-dDl+0.5*dDDld);
    dDl=-dDl;
    DDld=-DDld;
end
NDl0=norm(Dl);
Dl0=Dl;
if ((dDl)^2<=0.99*deltah*deltah*NDl0*NDl0)&&( NDl0>=0.1*lxd/deltah)
    %%
    %s在d与Dl张成的空间内
    %s=C1*d+C2*Dl;
    C2=deltah/(sqrt(NDl0*NDl0-(dDl/deltah)^2)   );
    C1=-C2*dDl/deltah/deltah;
    s=C1*d+C2*Dl;
    Dlxd=g+G*(xopt+d-x0);%Dlt(xopt+d0)
else
    %%
    %s在d0与Dl(xopt+d0)张成的空间内
    %s=C1*d+C2*Dl;
    Dlxd=g+G*(xopt+d-x0);%Dlt(xopt+d0);
    NDlxd=norm(Dlxd);
    dDlxd=d'*Dlxd;
    C2=deltah/(sqrt(NDlxd*NDlxd-(dDlxd/deltah)^2));
    C1=-C2*dDlxd/deltah/deltah;
    s=C1*d+C2*Dlxd;
end
theta=zeros(50,1);
for i=1:50
    theta(i)=2*(i-1)*pi/50;
end
%l(x+d)=cos*d'*Dl+sin*s'*Dl+0.5*cos^2*d'*DDl*d+0.5*sin^2*s'*DDl*s+sin*cos*s'*DDl*d
sDl=s'*Dl;
sDDld=s'*DDld;
DDls=DDlu(s);
sDDls=s'*DDls;
L=zeros(50,1);
for i=1:50
    C=cos(theta(i));
    S=sin(theta(i));
    L(i)=abs(C*dDl+S*sDl+0.5*(C*C*dDDld+S*S*sDDls)+S*C*sDDld);
end
lxdold=lxd;
[lxdnew,k]=max(L);
d=cos(theta(k))*d+sin(theta(k))*s;
DDld=DDlu(d);
times=1;
while (times<=N)&&(lxdnew>1.1*lxdold)
    times=times+1;
    lxdold=lxdnew;
    Dlxd=(1-cos(theta(k)))*Dl0+cos(theta(k))*Dlxd+sin(theta(k))*DDls;
    NDlxd=norm(Dlxd);
    dDlxd=d'*Dlxd;
    C2=deltah/(sqrt(NDlxd*NDlxd-(dDlxd/deltah)^2)   );
    C1=-C2*dDlxd/deltah/deltah;
    s=C1*d+C2*Dlxd;%新的s
    %l(x+d)=cos*d'*Dl+sin*s'*Dl+0.5*cos^2*d'*DDl*d+0.5*sin^2*s'*DDl*s+sin*cos*s'*DDl*d
    sDl=s'*Dl0;
    dDl=d'*Dl;
    dDDld=d'*DDld;
    sDDld=s'*DDld;
    DDls=DDlu(s);
    sDDls=s'*DDls;
    L=zeros(50,1);
    for i=1:50
        C=cos(theta(i));
        S=sin(theta(i));
        L(i)=abs(C*dDl+S*sDl+0.5*(C*C*dDDld+S*S*sDDls)+S*C*sDDld);
    end
    [lxdnew,k]=max(L);
    d=cos(theta(k))*d+sin(theta(k))*s;
    DDld=DDlu(d);
end
NORMD=deltah;
%%
xnew=xopt+d;
alpha=H(t,t);
w=zeros(m+n+1,1);
xx0=xnew-x0;
for i=1:m
    w(i)=0.5*((W(i,m+2:m+n+1)*xx0)^2);
end
w(m+1)=1;
w(m+2:m+n+1)=xx0;
Hw=H*w;
beta=0.5*(norm(xx0)^4)-w'*Hw;
tau=Hw(t);
tau2=tau*tau;
if abs(alpha*beta+tau2)<=0.8*tau2
    %%
    %BIGDEM
    Weitht=zeros(m-1,1);
    Weitht(t)=((XXopt(:,t)'*d)^2)/(NXX(t)^2*NORMD*NORMD);
    v=W(:,opt);
    if Weitht(t)<=0.99
        u=XXopt(:,t);
        NU=NXX(t);
    else
        for i=1:t-1
            Weitht(i)=((XXopt(:,i)'*d)^2)/(NXX(i)^2*NORMD*NORMD);
        end
        for i=t+1:m
            Weitht(i)=((XXopt(:,i)'*d)^2)/(NXX(i)^2*NORMD*NORMD);
        end
        Weitht(opt)=NaN;
        [~,k]=min(Weitht);
        u=XXopt(:,k);
        NU=NXX(k);
    end
    %s在d与u张成的空间内
    %s=C1*d+C2*u;
    ud=u'*d;
    C2=deltah/(sqrt(NU*NU-(ud/deltah)^2)   );
    C1=-C2*ud/deltah/deltah;
    s=C1*d+C2*u;
    DEN=zeros(50,1);
    xx1=xopt-x0;
    N1=norm(xx1);
    w=zeros(m+n+1,1);
    w(m+1)=1;
    for j=1:50
        dd=cos(theta(j))*d+sin(theta(j))*s;
        xnew=xopt+dd;
        xx0=xnew-x0;
        for i=1:m
            w(i)=0.5*((W(i,m+2:m+n+1)*xx0)^2);
        end
        w(m+2:m+n+1)=xx0;
        N0=norm(xx0);
        wv=w-v;
        Hwv=H*wv;
        DEN(j)= alpha*(0.5*N0^4-(xx1'*xx0)^2+0.5*N1^4     )  -alpha*wv'*Hwv      +(Hwv(t))^2     ;
    end
    [~,I]=max(abs(DEN));
    if I==1
        i1=50;i2=1;i3=2;
    else
        if I==50
            i1=49;i2=50;i3=1;
        else
            i1=I-1;i2=I;i3=I+1;
        end
    end
    q1=DEN(i1);q2=DEN(i2);q3=DEN(i3);
    theta1=2*(I-1)*pi/50;
    theta2=2*(I)*pi/50;
    theta3=2*(I+1)*pi/50;
    MI=[theta1*theta1,theta1,1;
        theta2*theta2,theta2,1;
        theta3*theta3,theta3,1];
    ce=MI\[q1;q2;q3];
    if abs(ce(1))>0
        mp=-ce(2)/(2*ce(1));
        if mp>theta1 && mp<theta3
            qmp=ce(1)*mp*mp+ce(2)*mp+ce(3);
            value=[q1;q2;q3;qmp];
            thetaC=[theta1;theta2;theta3;qmp];
            [~,qk]=max(abs(value));
            thetaq=thetaC(qk);
        else
            thetaq=theta2;
        end
    else
        thetaq=theta2;
    end
    d=cos(thetaq)*d+sin(thetaq)*s;%更新后的d
    xnew=xopt+d;
    w=zeros(m+n+1,1);
    xx0=xnew-x0;
    %xx1=xnew-x0;
    %N1=norm(xx1);
    for i=1:m
        w(i)=0.5*((W(i,m+2:m+n+1)*xx0)^2);
    end
    w(m+1)=1;
    w(m+2:m+n+1)=xx0;
    Hw=H*w;
    tau=Hw(t);
    %         N0=norm(xx0);
    wv=w-v;
    %         eta1=zeros(m,1);
    %         eta2=zeros(n,1);
    Hwv=H*wv;
    eta1=Hwv(1:m);
    eta2=Hwv(m+2:m+n+1);
    DDEN2=zeros(n,1);
    DDEN3=zeros(n,1);
    for i=1:m
        DDEN2=DDEN2+((tau*H(t,i)-alpha*eta1(i))*(xx0'*X(:,i)))*X(:,i);
    end
    DDEN2=2*DDEN2;
    for i=1:n
        DDEN3(i)=2*(tau*H(t,i+m+1)-alpha*eta2(i));
    end
    DDEN=2*alpha*(N1*N1*d+d'*xx1*xx0)+(DDEN2)+(DDEN3);
    NDDEN=norm(DDEN);
    times=1;
    while (times<=N)&&(((d'*DDEN)^2)<((1-10^-8)*deltah*NDDEN*NDDEN))
        times=times+1;
        u=DDEN;
        ud=u'*d;
        NU=NDDEN;
        C2=deltah/(sqrt(NU*NU-(ud/deltah)^2)   );
        C1=-C2*ud/deltah/deltah;
        s=C1*d+C2*u;
        w=zeros(m+n+1,1);
        w(m+1)=1;
        for j=1:50
            dd=cos(theta(j))*d+sin(theta(j))*s;
            xnew=xopt+dd;
            xx0=xnew-x0;
            for i=1:m
                w(i)=0.5*((W(i,m+2:m+n+1)*xx0)^2);
            end
            w(m+2:m+n+1)=xx0;
            N0=norm(xx0);
            wv=w-v;
            Hwv=H*wv;
            DEN(j)= alpha*(0.5*N0^4-(xx1'*xx0)^2+0.5*N1^4     )  -alpha*wv'*Hwv      +(Hwv(t))^2     ;
        end
        [~,I]=max(abs(DEN));
        if I==1
            i1=50;i2=1;i3=2;
        else
            if I==50
                i1=49;i2=50;i3=1;
            else
                i1=I-1;i2=I;i3=I+1;
            end
        end
        q1=DEN(i1);q2=DEN(i2);q3=DEN(i3);
        theta1=2*(I-1)*pi/50;
        theta2=2*(I)*pi/50;
        theta3=2*(I+1)*pi/50;
        MI=[theta1*theta1,theta1,1;
            theta2*theta2,theta2,1;
            theta3*theta3,theta3,1];
        ce=MI\[q1;q2;q3];
        if abs(ce(1))>0
            mp=-ce(2)/(2*ce(1));
            if mp>theta1 && mp<theta3
                qmp=ce(1)*mp*mp+ce(2)*mp+ce(3);
                value=[q1;q2;q3;qmp];
                thetaC=[theta1;theta2;theta3;qmp];
                [~,qk]=max(abs(value));
                thetaq=thetaC(qk);
            else
                thetaq=theta2;
            end
        else
            thetaq=theta2;
        end
        d=cos(thetaq)*d+sin(thetaq)*s;%更新后的d
        xnew=xopt+d;
        w=zeros(m+n+1,1);
        xx0=xnew-x0;
        for i=1:m
            w(i)=0.5*((W(i,m+2:m+n+1)*xx0)^2);
        end
        w(m+1)=1;
        w(m+2:m+n+1)=xx0;
        Hw=H*w;
        tau=Hw(t);
        wv=w-v;
        Hwv=H*wv;
        eta1=Hwv(1:m);
        eta2=Hwv(m+2:m+n+1);
        DDEN2=zeros(n,1);
        DDEN3=zeros(n,1);
        for i=1:m
            DDEN2=DDEN2+((tau*H(t,i)-alpha*eta1(i))*(xx0'*X(:,i)))*X(:,i);
        end
        DDEN2=2*DDEN2;
        for i=1:n
            DDEN3(i)=2*(tau*H(t,i+m+1)-alpha*eta2(i));
        end
        DDEN=2*alpha*(N1*N1*d+d'*xx1*xx0)+(DDEN2)+(DDEN3);
        NDDEN=norm(DDEN);
    end
end
%Xn(:,MOVE)=Xn(:,opt)+d;
xopt=Xn(:,opt);
xnew=xopt+d;%新生成的点
Fnew=F(xnew);
F_times=F_times+1;
Krho=Krho+1;
RATIO=1;
NEWUOAStep5();%To Setp5
end