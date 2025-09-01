a = [3.25 1.90 1]; %relative volatility
xf = [0.3 0.25 0.45];
xa_d = 0.98;
xc_d = 5*10^(-11);
xb_d = 0.02;
r = 2.7; %L/D = r = reflux ratio
F = 1;
xD = [0.98 0.02 5*10^(-11)];

%% part c(i)
xa_b =  0.005;

%D/B = k
k = (xf(1)-xa_b)/(xa_d-xf(1));
xb_b = xf(2) - k*(xb_d-xf(2));
xc_b = xf(3) - k*(xc_d-xf(3));

%V/D = q
q = r+1;

%V/B = s
s = k*q;

N=50;
%topdown
yr = zeros(N+1,length(xD));
xr = zeros(N,length(xD));

yr(1,:)=xD;
for i = 1:N
    xr(i,:)=yr(i,:)./a/(sum(yr(i,:)./a)); %VLE
    yr(i+1,:)= r/(r+1)*xr(i,:)+xD/(r+1); %CMB
end

%bottomup
xB = [0.005 0.350 0.645];

ys = zeros(N,length(xB));
xs = zeros(N+1,length(xB));

xs(1,:)=xB;
for j = 1:N
    ys(j,:)=xs(j,:).*a/(sum(xs(j,:).*a)); %VLE
    xs(j+1,:)= s/(s+1)*ys(j,:)+xB/(s+1); %CMB 
end

hold on
plot(xr(:,1),xr(:,2),'b*')
plot(xs(:,1),xs(:,2),'k*')
plot([xB(1,1),xD(1,1)],[xB(1,2),xD(1,2)],'r')
xlabel('xa'); ylabel('xb'); title('c(i)',FontWeight='normal');
%legend('xr','xs','line');
hold off
