clc
clear all

%% Constants:

% Antoine Constants [B T X]
A = [6.905 6.954 6.99];
B = [1211.03 1344.8 1453.43];
C = [220.79 219.48 215.31];
P_tot = 1.1*750.061; %pressure total
F = 100; % kmol/hr feed flow rate
zf = [1,1,1]/3; % initial feed composition
lhv_bz = 33900; % kJ/kmol
lhv_tl = 35000; % kJ/kmol
rho_L = 786.3; % kg/m³
rho_V = 2.577; % kg/m³
% 1 = benzene; 2 = toluene; 3 = xylene;

%% Relative volatility cal

% Benzene & Xylene 
x_grid = linspace(0,1,100);
for i = 1:length(x_grid)
    fun1 = @(T) (P_tot-x_grid(i)*10^(A(1)-B(1)/(C(1)+T))-(1-x_grid(i))*10^(A(3)-B(3)/(C(3)+T)));
    T_bx(i) = fzero(fun1,30);
    y_bx(i) = x_grid(i)*10^(A(1)-B(1)/(C(1)+T_bx(i)))/P_tot;
end

% Least square curve fit 
fun3 = @(alpha,x) (alpha.*x./(alpha.*x+(1-x))); 
alpha(1) = lsqcurvefit(fun3,1,x_grid,y_bx);

% Toluene & Xylene
for i = 1:length(x_grid)
    fun2 = @(T) (P_tot-x_grid(i)*10^(A(2)-B(2)/(C(2)+T))-(1-x_grid(i))*10^(A(3)-B(3)/(C(3)+T)));
    T_tx(i) = fzero(fun2,30);
    y_tx(i) = x_grid(i)*10^(A(2)-B(2)/(C(2)+T_tx(i)))/P_tot;
end

alpha(2) = lsqcurvefit(fun3,1,x_grid,y_tx);
alpha(3) = 1;
global a
a = alpha;

plot(x_grid,y_bx,x_grid,y_tx,x_grid,x_grid)
legend('Benzene-Xylene','Toluene-Xylene','y=x',Location='southeast');
xlabel('x (liq mole frac)'); ylabel('y (vap mole frac)');
axis([0 1 0 1])

%% Feasibility limit

a = alpha;  
x_d1_grid = 1:-0.01:zf(1);
y2_sol = zeros(length(x_d1_grid),1);
tol = 1e-1;

for i = 2:length(x_d1_grid)
    f = @(y2) str_locus(a,zf,x_d1_grid(i),y2);
    y2_guess = 1-abs(x_d1_grid);
    for j = 1:length(y2_guess)
        sol_return(j) = fzero(f,y2_guess(j));
        k = abs(sol_return(j)-y2_sol(i-1));
        if(k<tol)
            y2_sol(i) = sol_return(j);
        break
        end
    end
end

y2_sol=y2_sol';

% Maximum impurity
x_d1 = 0.99;
[~, idx] = min(abs(x_d1_grid - x_d1));
feasibility_limit = y2_sol(idx);

% Rmin of heavy key
x_d2_grid = linspace(1e-8,feasibility_limit,51);
x_d3_grid = 1-x_d2_grid-x_d1;
x_b1 = 0.01;
k1 = (zf(1)-x_b1)/(x_d1-zf(1)); % lever rule

%B = bottom flow; D = distillate flow;
B = F/(k1+1);
D = k1*B;
R_min = zeros(1, length(x_d2_grid));  

for j = 1:length(x_d2_grid)
    xD = [x_d1, x_d2_grid(j), x_d3_grid(j)];
    xB = [x_b1, zf(2) - k1 * (xD(2) - zf(2)), zf(3) - k1 * (xD(3) - zf(3))];
    b = linspace(0.95, 0.01, length(x_d2_grid));
    area = zeros(length(b),1);
    xpr_guess = [0, 1];
    xps_guess = [1, 0];

    for p = 1:length(b)
        r = b(p) / (1 - b(p));
        s = (r + 1) * k1;
        xD_trim = xD(1:2);
        xB_trim = xB(1:2);
        a_trim = a(1:2);
        f = @(x) rctf_sadle(xD_trim, r, a_trim, x);
        f2 = @(x) str_node(xB_trim, s, a_trim, x);
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

%disp('Final R_min values:');
%disp(R_min);
R_operation = 1.25 .* R_min;

for j = 1:length(R_operation)
    R = R_operation(j);
    N = 50;
    xD = [x_d1, x_d2_grid(j), x_d3_grid(j)];
    xB = [x_b1, zf(2) - k1 * (xD(2) - zf(2)), zf(3) - k1 * (xD(3) - zf(3))];
    xa_d = xD(1);
    xa_b = xB(1);
    k1_val = (zf(1) - xa_b) / (xa_d - zf(1)); % Use scalar for k1 calculation

    B = F / (k1_val + 1);
    D = k1_val * B;

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
        idx = j + 1;
    end
    idx_s = 0; 
    
    for j = 1:N
        diff = -(xs(j,1) - zf(1));
        if diff < 0
           idx_s = j;
           break
        end
        idx_s = j + 1;
    end

    N_rctf(j) = idx;
    N_str(j) = idx_s;
    N_total(j) = N_str(j) + N_rctf(j);
end

[~, idx_opt] = min(N_total);
xD2_opt = x_d2_grid(idx_opt);  
R_opt = R_min(idx_opt);
N_total