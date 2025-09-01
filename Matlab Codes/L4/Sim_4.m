clc
clear all
%% Relative volatility cal

% Antoine Constants [B T X]
A = [6.905 6.954 6.99];
B = [1211.03 1344.8 1453.43];
C = [220.79 219.48 215.31];

P=1.1*750.061;

% Benzene and Xylene 
x_grid = linspace(0,1,100);
for i = 1:length(x_grid)
    fun1 = @(T) (P-x_grid(i)*10^(A(1)-B(1)/(C(1)+T))-(1-x_grid(i))*10^(A(3)-B(3)/(C(3)+T)));
    T_bx(i)=fzero(fun1,30);
    y_bx(i)=x_grid(i)*10^(A(1)-B(1)/(C(1)+T_bx(i)))/P;
end

%non-linear least square curve fit 
fun1= @(alpha,x) (alpha.*x./(alpha.*x+(1-x))); 
alpha(1) = lsqcurvefit(fun1,1,x_grid,y_bx);


% Toluene and Xylene
for i=1:length(x_grid)
    fun1 = @(T)(P-x_grid(i)*10^(A(2)-B(2)/(C(2)+T))-(1-x_grid(i))*10^(A(3)-B(3)/(C(3)+T)));
    T_tx(i)=fzero(fun1,30);
    y_bx(i)=x_grid(i)*10^(A(2)-B(2)/(C(2)+T_tx(i)))/P;
end

fun1=@(alpha,x)(alpha.*x./(alpha.*x+(1-x)));%non-linear least square curve fit
alpha(2)=lsqcurvefit(fun1,1,x_grid,y_bx);
alpha(3)=1;
global a
a=alpha;

plot(x_grid,y_bx)

%% Feasibility limit

F=100;
zf=[1,1,1]/3; % feed composition
LHV_BENZ = 33900; %kJ/kmol
LHV_TOL = 35000; %kJ/kmol
rho_L = 786.3; %kg/m³
rho_V = 2.577; %kg/m³
a=alpha;  
xD1_array=1:-0.01:zf(1);
y2_sol=zeros(length(xD1_array),1);
tol=1e-1;
for i=2:length(xD1_array)
    f=@(y2) stripping_locus(a,zf,xD1_array(i),y2);
    y2_guess=1-abs(xD1_array);
    for j=1:length(y2_guess)
        sol_return(j)=fsolve(f,y2_guess(j));
        k=abs(sol_return(j)-y2_sol(i-1));
        if(k<tol)
            y2_sol(i)=sol_return(j);
        break
        end
    end
end
y2_sol=transpose(y2_sol);

%need to find maximum impurity corresponding to given xD 
xD1 = 0.99;
[~, idx] = min(abs(xD1_array - xD1));
feasibility_limit = y2_sol(idx);

% Rmin of heavy key
xD2_array=linspace(1e-8,feasibility_limit,51);
xD3_array=1-xD2_array-xD1;
xB1=0.01;
m=(zf(1)-xB1)/(xD1-zf(1));%LEVER RULE 

% GETTING DISTILLATE AND BOTTOM FLOWS 
  B=F/(m+1);
  D=m*B;
R_min = zeros(1, length(xD2_array));  
for j = 1:length(xD2_array)
    xD = [xD1, xD2_array(j), xD3_array(j)];
    xB = [xB1, zf(2) - m * (xD(2) - zf(2)), zf(3) - m * (xD(3) - zf(3))];

    b = linspace(0.95, 0.01, length(xD2_array));
    area = zeros(length(b), 1);

    xpr_guess = [0, 1];
    xps_guess = [1, 0];

    for p = 1:length(b)
        r = b(p) / (1 - b(p));
        s = (r + 1) * m;
        
        xD_trim = xD(1:2);
        xB_trim = xB(1:2);
        a_trim = a(1:2);

        f = @(x) rectifying_sadle(xD_trim, r, a_trim, x);
        f2 = @(x) stripping_node(xB_trim, s, a_trim, x);

        xp_r = fsolve(f, xpr_guess);
        xp_s = fsolve(f2, xps_guess);

        xpr_guess = xp_r;
        xps_guess = xp_s;

        e1 = (zf(1:2) - xp_r)';
        e2 = (xp_s - xp_r)';
        area(p) = det([e1 e2]);
    end

    idx = find(area < 0, 1);
  
    disp(['Iteration j = ', num2str(j)]);
    disp(['Index found: ', num2str(idx)]);
    
    if ~isempty(idx)
        b_req = b(idx);
        R_min(j) = b_req / (1 - b_req);
    end
end

disp('Final R_min values:');
disp(R_min);
R_operation=1.25.*R_min;
for j=1:length(R_operation)
    R=R_operation(j);
    N = 50;
 xD = [xD1, xD2_array(j), xD3_array(j)];
 xB = [xB1, zf(2) - m * (xD(2) - zf(2)), zf(3) - m * (xD(3) - zf(3))];
m = @(xa_b) (zf(1) - xa_b) / (xa_d - zf(1));
xa_b =xB(1);

B = F / (m(xa_b) + 1);
D = m(xa_b) * B;

V = (R + 1) * D;
L_rectifying = R * D;
S = V / B;
L_stripping = (S + 1) * B;

yr = zeros(N+1, 3);
xr = zeros(N, 3);
xs = zeros(N+1, 3);
ys = zeros(N, 3);
yr(1, :) = xD;

for i = 1:N
    xr(i, :) = yr(i, :) ./ a / sum(yr(i, :) ./ a);
    yr(i+1, :) = (R / (R+1)) * xr(i, :) + xD / (R+1);
end

xs(1, :) = xB;
for i = 1:N
    ys(i, :) = xs(i, :) .* a / sum(xs(i, :) .* a);
    xs(i+1, :) = (S / (S+1)) * ys(i, :) + xB / (S+1);
end 

tol = 1e-4;

idx = 0; 

for j = 1:N
    diff = (xr(j,1) - zf(1));
    if diff < 0
       idx = j;
       break
    end
    idx=j+1;
end
idx_s = 0; 

for j =1:N
    diff = -(xs(j,1) - zf(1));
    if diff < 0
       idx_s = j;
       break
    end
    idx_s=j+1;
end

    N_rectifying(j) = idx;
    N_stripping(j) = idx_s;
    N_total(j)=N_stripping(j)+N_rectifying(j);
end
[~,idx_opt]=min(N_total);
xD2_opt=xD2_array(idx_opt); % optimum heavy key impurity to operate at 
R_opt=R_min(idx_opt);

% PART B 
R_array=R_opt:0.1:10*R_opt;
for j=1:length(R_array)
       R=R_array(j);
    N = 50;
 xD = [xD1, xD2_array(j), xD3_array(j)];
 xB = [xB1, zf(2) - m * (xD(2) - zf(2)), zf(3) - m * (xD(3) - zf(3))];
m = @(xa_b) (zf(1) - xa_b) / (xa_d - zf(1));
xa_b =xB(1);

B = F / (m(xa_b) + 1);
D = m(xa_b) * B;

V = (R + 1) * D;
L_rectifying = R * D;
S = V / B;
L_stripping = (S + 1) * B;

yr = zeros(N+1, 3);
xr = zeros(N, 3);
xs = zeros(N+1, 3);
ys = zeros(N, 3);
yr(1, :) = xD;

for i = 1:N
    xr(i, :) = yr(i, :) ./ a / sum(yr(i, :) ./ a);
    yr(i+1, :) = (R / (R+1)) * xr(i, :) + xD / (R+1);
end

xs(1, :) = xB;
for i = 1:N
    ys(i, :) = xs(i, :) .* a / sum(xs(i, :) .* a);
    xs(i+1, :) = (S / (S+1)) * ys(i, :) + xB / (S+1);
end 

tol = 1e-4;

idx = 0; 

for j = 1:N
    diff = (xr(j,1) - zf(1));
    if diff < 0
       idx = j;
       break
    end
    idx=j+1;
end
idx_s = 0; 

for j =1:N
    diff = -(xs(j,1) - zf(1));
    if diff < 0
       idx_s = j;
       break
    end
    idx_s=j+1;
end

    N_rect(j) = idx;
    N_strip(j) = idx_s;
    N_tot(j)=N_strip(j)+N_rect(j);
    % COST CALCULATIONS 
    Q_reb =((V-D)*LHV_TOL)/3600;
    Q_cond =(V*LHV_BENZ)/3600;
    FLV=(L_rectifying/V)*(rho_V/rho_L)^0.5;
end
