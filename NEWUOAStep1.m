function [] = NEWUOAStep1()
%NEWUOAStep1 对应NEWUOA算法的第一步
%   在第一步中，会根据初始点生成初始的插值节点，以及计算它们处的函数值，同时计算此时的最小值的Fopt以及最小值点xopt
%   Xn是插值点矩阵，W是插值问题的矩阵，H是W的逆（不用求逆算法得）,Fn是插值处目标函数的函数值
%   最初始的插值函数Q(x)=c+g'*(x-x0)+0.5*(x-x0)'*G*(x-x0);
%   W=[A,X';X,0]; H=[s.*Z*Z';Theta';Theta,Upsilon];
global Xn Fn n m F_times rho_beg x0 opt xb c g Gamma gamma H F rho delta Krho D3 QF3
gamma=zeros(m,1);
Xn=zeros(n,m);
Fn=zeros(m,1);
x0=xb;
Xn(:,1)=x0;
Fn(1)=F(Xn(:,1));
I=eye(n);
Gamma=zeros(n,n);
c=Fn(1);
g=zeros(n,1);
%%
%设置参数
rho=rho_beg;
delta=rho;
Krho=0;
D3=zeros(3,1);
QF3=zeros(3,1);
%%
%生成插值点
if m>=2*n+1
    Delta=zeros(n,1);
    for i=1:n
        Xn(:,i+1)=x0+rho*I(:,i);
        Fn(i+1)=F(Xn(:,i+1));
        Xn(:,i+n+1)=x0-rho*I(:,i);
        Fn(i+n+1)=F(Xn(:,i+n+1));
        if Fn(i+n+1)>=Fn(i+1)
            Delta(i)=1;
        else
            Delta(i)=-1;
        end
        g(i)=(Fn(i+1)-Fn(i+n+1))/(2*rho);
        Gamma(i,i)=(Fn(i+1)+Fn(i+n+1)-2*Fn(1))/(rho*rho);
    end
    
    for i=2*n+2:1:m
        j=floor((i-n-2)/n);
        p=i-n-1-j*n;
        if p+j>n
            q=p+j-n;
        else
            q=p+j;
        end
        Xn(:,i)=x0+rho*Delta(p)*I(:,p)+rho*Delta(q)*I(:,q);
        Fn(i)=F(Xn(:,i));
        Gamma(p,q)=( Fn(1)-Fn(1+p+(1-Delta(p))*n/2)-Fn(1+q+(1-Delta(q))*n/2)+Fn(i) )/(Delta(p)*Delta(q)*rho*rho);
        Gamma(q,p)=Gamma(p,q);
    end
else
    for i=1:n
        Xn(:,i+1)=x0+rho*I(:,i);
        Fn(i+1)=F(Xn(:,i+1));
    end
    for i=n+1:m-1
        Xn(:,i+1)=x0-rho*I(:,i-n);
        Fn(i+1)=F(Xn(:,i+1));
        g(i-n)=(Fn(i-n+1)-Fn(i+1))/(2*rho);
        Gamma(i-n,i-n)=(Fn(i-n+1)+Fn(i+1)-2*Fn(1))/(rho*rho);
    end
    for i=m:2*n
        g(i-n)=(Fn(i-n+1)-Fn(1))/(rho);
    end
end
F_times=m;
%%
%计算W
% XX0=Xn-x0;
% A=((XX0'*XX0).^2)/2;
% X=[ones(1,m);XX0];
% W=zeros(m+n+1,m+n+1);
% W(1:m,1:m)=A;
% W(m+1:m+n+1,1:m)=X;
% W(1:m,m+1:m+n+1)=X';
%%
%计算W的逆H
Theta=zeros(n+1,m);
Theta(1,1)=1;
Upsilon=zeros(n+1,n+1);
if m>=2*n+1
    for i=2:n+1
        Theta(i,i)=1/(2*rho);
        Theta(i,i+n)=-1/(2*rho);
    end
else
    for i=2:m-n
        Theta(i,i)=1/(2*rho);
        Theta(i,i+n)=-1/(2*rho);
    end
    for i=m-n+1:n+1
        Theta(i,1)=-1/rho;
        Theta(i,i)=1/rho;
        Upsilon(i,i)=-0.5*rho*rho;
    end
end
Z=zeros(m,m-n-1);
s=ones(1,m-n-1);
if m<=2*n+1
    for k=1:m-n-1
        Z(1,k)=-sqrt(2)/(rho*rho);
        Z(k+1,k)=-Z(1,k)/2;
        Z(k+n+1,k)=Z(k+1,k);
    end
else
    for k=1:n
        Z(1,k)=-sqrt(2)/(rho*rho);
        Z(k+1,k)=-Z(1,k)/2;
        Z(k+n+1,k)=Z(k+1,k);
    end
    for k=n+1:m-n-1
        i=k+n+1;
        j=floor((i-n-2)/n);
        p=i-n-1-j*n;
        if p+j>n
            q=p+j-n;
        else
            q=p+j;
        end
        if Delta(p)==1
            pt=p+1;
        else
            pt=p+1+n;
        end
        if Delta(q)==1
            qt=q+1;
        else
            qt=q+1+n;
        end
        Z(1,k)=1/(rho*rho);
        Z(pt,k)=-Z(1,k);
        Z(qt,k)=Z(pt,k);
        Z(k+n+1,k)=Z(1,k);
    end
end
H=zeros(m+n+1,m+n+1);
H(1:m,1:m)=Z*Z';
H(1:m,m+1:m+n+1)=Theta';
H(m+1:m+n+1,1:m)=Theta;
H(m+1:m+n+1,m+1:m+n+1)=Upsilon;
%%
%计算二次模型系数c,g,lambda
% R=zeros(m+n+1,1);
% R(1:m)=Fn;
% Lcg=H*R;
% lambda=Lcg(1:m);
% c=Lcg(m+1);
% g=Lcg(m+2:m+n+1);
%%
[~,opt]=min(Fn);
%%
NEWUOAStep2();
end

