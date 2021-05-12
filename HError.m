function [HE,W,Hreal] = HError()
%HERROR æÿ’ÛHµƒŒÛ≤Ó
%   ≤‚ ‘”√∫Ø ˝

global H Xn m n x0 
X=zeros(n+1,m);
X(1,:)=ones(1,m);
X(2:n+1,:)=Xn-x0;
A=zeros(m,m);
for i=1:m
    for j=1:m
        A(i,j)=0.5*((Xn(:,i)-x0)'*(Xn(:,j)-x0))^2;
    end
end
W=zeros(m+n+1,m+n+1);
W(1:m,1:m)=A;
W(1:m,m+1:m+n+1)=X';
W(m+1:m+n+1,1:m)=X;
Hreal=inv(W);
HE=sum(sum(abs(H-Hreal)));
end

