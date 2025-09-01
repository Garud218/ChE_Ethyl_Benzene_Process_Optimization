function F = eqn1(x)
global x1 z a

y1p = a(1)*x1/(a(1)*x1+a(2)*x+a(3)*(1-x1-x));
y2p = a(2)*x/(a(1)*x1+a(2)*x+a(3)*(1-x1-x));

F = (y2p-x)/(y1p-x1) - (z(2)-x)/(z(1)-x1);
end