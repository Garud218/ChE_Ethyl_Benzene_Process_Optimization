global x1 z a
z = [0.3 0.3 0.4];
a = [12.67 5.35 1];

xvar = linspace(0.99,0.01,100); x2 = zeros(1,100);
for i = 1:100
    x1 = xvar(i);
    xguess = 1-x1;
    x2(i) = fsolve(@eqn1,xguess);
end

figure
plot(xvar(1:71),x2(1:71))
hold on
plot(0:0.01:0.5,0:0.01:0.5)
plot(0:0.01:1,1:-0.01:0)
plot(0.3,0.3,'k*')
hold off
