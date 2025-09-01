a = [4 2 1]; %relative volatility
xf = [1/3 1/3 1/3];
xa_d = 0.99;
xc_d = 1e-10;
xb_d = 1-xa_d-xc_d;
r = 2.1566*1.25; %L/D = r = reflux ratio
F = 1; r0=3;
xD = [xa_d xb_d xc_d];

xa_b =  0.01;

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
xB = [xa_b xb_b xc_b];

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
xlabel('xa'); ylabel('xb'); title('1.25*R_{min} plot',FontWeight='normal');
%legend('xr','xs','line');
hold off 

%b=0.6832
%r=2.1566
b= [0.95:-0.01:0.1]';
M = length(b);
area = zeros(M,1);
r1 = area;
xpr = [0 1 0]';
xps = [1 0 0]';
s0=(r0+1)*k;
for i=1:M
    r1 = b(i)/(1-b(i));
    s1 = (r1+1)*k;
    xpr = crf(xD,r1,xpr);
    xps = csf(xB,s1,xps);
    e1 = xf(1:2)-xpr(1:2); 
    e2 = xps(1:2)-xpr(1:2);
    area(i) = det([e1, e2]);
end

figure
plot(b,area)