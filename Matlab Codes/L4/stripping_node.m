function F=stripping_node(xB,s,a,x)
k1=s/(s+1);
k2=k1/s;
F(1)=x(1)-k1*(a(1)*x(1)/(1+(a(1)-1)*x(1)+(a(2)-1)*x(2))-k2*xB(1));
F(2)=x(2)-k1*(a(2)*x(2)/(1+(a(1)-1)*x(1)+(a(2)-1)*x(2)))-k2*xB(2);
end

