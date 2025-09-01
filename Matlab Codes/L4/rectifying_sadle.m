function F=rectifying_sadle(xD,r,a,x)
k1=r/(r+1);
k2=k1/r;
F(1)=a(1)*x(1)/(1+(a(1)-1)*x(1)+(a(2)-1)*x(2))-k1*x(1)-k2*(xD(1));
F(2)=a(2)*x(2)/(1+(a(1)-1)*x(1)+(a(2)-1)*x(2))-k1*x(2)-k2*(xD(2));
end

