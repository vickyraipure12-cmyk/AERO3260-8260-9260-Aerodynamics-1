%% AERO3260 Assignment 1 - Task 3: Thin Aerofoil Theory
% SID = 550361509
% N7 = 5, N8 = 0, N9 = 9
% Camber line + analytical Fourier coefficients + 10-panel LVM
% Final version for report calculations.

clear; clc; close all;

%% -------------------- INPUTS --------------------
c = 1.0;                    % chord [m]
Uinf = 1.0;                 % freestream velocity
N7 = 5; N8 = 0; N9 = 9;
alpha_deg = N7;
alpha = deg2rad(alpha_deg);

fprintf('================ TASK 3: THIN AEROFOIL THEORY ================\n');
fprintf('SID = 550361509\n');
fprintf('N7 = %d, N8 = %d, N9 = %d\n',N7,N8,N9);
fprintf('c = %.4f m, Uinf = %.4f, alpha = %.4f deg\n\n',c,Uinf,alpha_deg);

%% -------------------- CAMBER LINE --------------------
% Given by the assignment:
% y/c = 3(N8+N9)(c-x)^2*x/(200*c^3)
% For N8+N9 = 9 and c = 1:
% y = 0.135*x*(1-x)^2

K = 3*(N8+N9)/200;
yc_fun = @(x) K*x.*(c-x).^2/c^3;
dyc_dx_fun = @(x) K*(c-x).*(c-3*x)/c^3;

x_plot = linspace(0,c,2001);
y_plot = yc_fun(x_plot);

% Maximum camber (analytical location for this polynomial)
xmax = c/3;
ymax = yc_fun(xmax);

fprintf('Camber line: y/c = 3(N8+N9)(c-x)^2*x/(200c^3)\n');
fprintf('For this SID: y/c = 0.135 (x/c)(1-x/c)^2\n');
fprintf('Maximum camber: x/c = %.8f, y/c = %.8f\n\n',xmax/c,ymax/c);

figure('Color','w','Position',[100 100 1000 650]);
plot(x_plot/c,y_plot/c,'k-','LineWidth',2.2); hold on;
plot(xmax/c,ymax/c,'ko','MarkerFaceColor','k','MarkerSize',6);
grid on; box on;
axis equal;
xlim([0 1]);
ylim([-0.02 0.08]);
xlabel('x/c','FontSize',13);
ylabel('y/c','FontSize',13);
title('Task 3(a): Aerofoil Camber Line','FontSize',14,'FontWeight','bold');
legend('Camber line','Maximum camber','Location','best');
set(gca,'FontSize',12,'LineWidth',1);

%% -------------------- ANALYTICAL FOURIER COEFFICIENTS --------------------
% Thin-aerofoil substitution:
% x/c = (1-cos(theta))/2, theta in [0,pi]

x_theta = @(th) c*(1-cos(th))/2;
slope_theta = @(th) dyc_dx_fun(x_theta(th));

I0 = integral(@(th) slope_theta(th).*(cos(th)-1),0,pi,'AbsTol',1e-12,'RelTol',1e-12);
A0_integral_term = I0/pi;

A1 = (2/pi)*integral(@(th) slope_theta(th).*cos(th),0,pi,...
    'AbsTol',1e-12,'RelTol',1e-12);

A2 = (2/pi)*integral(@(th) slope_theta(th).*cos(2*th),0,pi,...
    'AbsTol',1e-12,'RelTol',1e-12);

% At zero lift: CL = pi(2A0+A1) = 0.
% A0 = alpha + (integral term as defined in the assignment convention).
% Therefore alpha_0 = I0/pi - A1/2.
alpha0 = A0_integral_term - A1/2;
alpha0_deg = rad2deg(alpha0);

% At arbitrary alpha, A0(alpha) = alpha + I0/pi
A0_at_alpha = alpha + A0_integral_term;
CL_analytical = pi*(2*A0_at_alpha + A1);
CL_from_zero_lift = 2*pi*(alpha-alpha0);

Cm_c4_analytical = -(pi/4)*(A1-A2);

fprintf('================ ANALYTICAL THIN-AEROFOIL RESULTS ================\n');
fprintf('Integral I0/pi = %.12f\n',A0_integral_term);
fprintf('A0 at alpha = %.8f rad\n',A0_at_alpha);
fprintf('A1 = %.12f\n',A1);
fprintf('A2 = %.12f\n',A2);
fprintf('Zero-lift angle alpha0 = %.12f rad = %.8f deg\n',alpha0,alpha0_deg);
fprintf('CL (from pi(2A0+A1)) = %.12f\n',CL_analytical);
fprintf('CL (from 2*pi*(alpha-alpha0)) = %.12f\n',CL_from_zero_lift);
fprintf('Cm_c/4 = %.12f\n\n',Cm_c4_analytical);

%% -------------------- 10-PANEL LUMPED VORTEX MODEL --------------------
% Assignment requirement:
% 10 equal-length panels in x; vortex at quarter point;
% boundary condition at three-quarter point.

Np = 10;
Dx = c/Np;

x_vortex = ((0:Np-1)+0.25)*Dx;
x_control = ((0:Np-1)+0.75)*Dx;

y_vortex = yc_fun(x_vortex);
y_control = yc_fun(x_control);

slope_control = dyc_dx_fun(x_control);
normal_mag = sqrt(1+slope_control.^2);
nx = -slope_control./normal_mag;
ny = 1./normal_mag;

% Influence matrix: point vortex convention
% u = -Gamma/(2*pi) * dy/r^2
% v =  Gamma/(2*pi) * dx/r^2
A = zeros(Np,Np);

for i = 1:Np
    for j = 1:Np
        dxij = x_control(i)-x_vortex(j);
        dyij = y_control(i)-y_vortex(j);
        r2 = dxij^2 + dyij^2;
        uij = -dyij/(2*pi*r2);
        vij =  dxij/(2*pi*r2);
        A(i,j) = uij*nx(i) + vij*ny(i);
    end
end

% No-penetration: (Vinf + Vinduced).n = 0
rhs = -(Uinf*cos(alpha)*nx + Uinf*sin(alpha)*ny);
Gamma_panel = A\rhs;

% With the above vortex sign convention, positive lift corresponds to
% -sum(Gamma_panel). The total circulation is therefore:
Gamma_total = sum(Gamma_panel);
CL_LVM = -2*Gamma_total/(Uinf*c);

% Moment about quarter chord from the discrete vortex distribution.
Cm_c4_LVM = (2/(Uinf*c^2))*sum(Gamma_panel.*(x_vortex-c/4));

% Residual of boundary condition
BC_residual = A*Gamma_panel-rhs;

% Percentage differences
CL_diff_pct = abs(CL_LVM-CL_analytical)/abs(CL_analytical)*100;
Cm_diff_pct = abs(Cm_c4_LVM-Cm_c4_analytical)/abs(Cm_c4_analytical)*100;

fprintf('================ 10-PANEL LVM RESULTS ================\n');
fprintf('Number of panels = %d\n',Np);
fprintf('Panel length = %.8f m\n',Dx);
fprintf('Vortex location = quarter point of each panel\n');
fprintf('Control location = three-quarter point of each panel\n');
fprintf('Total discrete circulation = %.12f\n',Gamma_total);
fprintf('CL_LVM = %.12f\n',CL_LVM);
fprintf('Cm_c/4,LVM = %.12f\n',Cm_c4_LVM);
fprintf('Maximum BC residual = %.3e\n',max(abs(BC_residual)));
fprintf('CL difference from analytical = %.4f %%\n',CL_diff_pct);
fprintf('Cm difference from analytical = %.4f %%\n\n',Cm_diff_pct);

%% -------------------- LVM CIRCULATION DISTRIBUTION --------------------
figure('Color','w','Position',[120 120 1000 650]);
bar(x_vortex/c,Gamma_panel,'FaceColor',[0.2 0.2 0.2],'EdgeColor','k');
grid on; box on;
xlabel('x/c','FontSize',13);
ylabel('\Gamma_j / (U_\infty c)','FontSize',13);
title('Task 3(d): 10-Panel Lumped Vortex Model Circulation','FontSize',14,'FontWeight','bold');
set(gca,'FontSize',12,'LineWidth',1);

%% -------------------- COMPARISON TABLE IN COMMAND WINDOW --------------------
fprintf('================ FINAL TASK 3 VALUES ================\n');
fprintf('Analytical alpha0      = %.8f deg\n',alpha0_deg);
fprintf('Analytical CL          = %.8f\n',CL_analytical);
fprintf('Analytical A1          = %.8f\n',A1);
fprintf('Analytical A2          = %.8f\n',A2);
fprintf('Analytical Cm_c/4      = %.8f\n',Cm_c4_analytical);
fprintf('10-panel LVM CL        = %.8f\n',CL_LVM);
fprintf('10-panel LVM Cm_c/4    = %.8f\n',Cm_c4_LVM);
fprintf('======================================================\n');
