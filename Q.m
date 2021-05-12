function [value] = Q(x)
%Q 插值函数
%测试用函数
global c g Gamma gamma x0 Xn m
G=Gamma;
for i=1:m
    G=G+gamma(i)*((Xn(:,i)-x0)*(Xn(:,i)-x0)');
end
value=c+(x-x0)'*g+0.5*(x-x0)'*G*(x-x0);

end

