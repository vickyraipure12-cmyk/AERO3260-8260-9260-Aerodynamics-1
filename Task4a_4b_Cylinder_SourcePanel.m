%% AERO3260 - Aerodynamics 1
% Assignment 1 - Task 4(a) and 4(b)
% FINAL: 64-panel constant-source panel method for a circular cylinder
%
% Task 4(a): Streamlines around a unit-radius cylinder
% Task 4(b): Cp comparison with Cp = 1 - 4*sin(theta)^2
%
% The velocity-field calculation below uses the same constant-source
% panel formulation as the surface solution, with the standard local
% source-panel velocity:
%
%   u_local = lambda/(4*pi) * ln(r1^2/r2^2)
%   v_local = lambda/(2*pi) * (theta2-theta1)
%
% For the selected CCW geometry and outward-normal convention, the
% self-induced source-panel velocity is u_local = 0, v_local = -1/2.

clear;
clc;
close all;

fprintf('\n============================================================\n');
fprintf(' AERO3260 - TASK 4(a) AND 4(b)\n');
fprintf(' 64-PANEL CONSTANT-SOURCE CYLINDER METHOD\n');
fprintf('============================================================\n\n');

%% 1. Parameters
U_inf = 1.0;
alpha = 0.0;
R = 1.0;
N = 64;

fprintf('Freestream velocity U_inf = %.4f\n',U_inf);
fprintf('Cylinder radius R         = %.4f\n',R);
fprintf('Number of panels N        = %d\n',N);
fprintf('Angle of attack           = %.2f deg\n\n',rad2deg(alpha));

%% 2. Cylinder geometry - counter-clockwise
theta_nodes = linspace(0,2*pi,N+1);

x_nodes = R*cos(theta_nodes);
y_nodes = R*sin(theta_nodes);

x1 = x_nodes(1:N);
y1 = y_nodes(1:N);
x2 = x_nodes(2:N+1);
y2 = y_nodes(2:N+1);

xc = 0.5*(x1+x2);
yc = 0.5*(y1+y2);

S = hypot(x2-x1,y2-y1);
beta = atan2(y2-y1,x2-x1);

% Outward normal for CCW geometry
nx = sin(beta);
ny = -cos(beta);

% Tangential direction
tx = cos(beta);
ty = sin(beta);

%% 3. Source-panel influence matrix
A = zeros(N,N);

for i = 1:N
    for j = 1:N

        dx = xc(i)-x1(j);
        dy = yc(i)-y1(j);

        % Local coordinates relative to panel j
        X =  dx*cos(beta(j)) + dy*sin(beta(j));
        Y = -dx*sin(beta(j)) + dy*cos(beta(j));

        if i == j
            % Self influence for selected convention
            uLocal = 0;
            vLocal = -0.5;
        else
            r1sq = X^2 + Y^2;
            r2sq = (X-S(j))^2 + Y^2;

            uLocal = (1/(4*pi))*log(r1sq/r2sq);

            theta1 = atan2(Y,X);
            theta2 = atan2(Y,X-S(j));

            vLocal = (1/(2*pi))*(theta2-theta1);
        end

        % Local -> global
        uInd = uLocal*cos(beta(j)) - vLocal*sin(beta(j));
        vInd = uLocal*sin(beta(j)) + vLocal*cos(beta(j));

        % Normal-velocity influence coefficient
        A(i,j) = uInd*nx(i) + vInd*ny(i);
    end
end

%% 4. Solve source strengths
u_inf = U_inf*cos(alpha);
v_inf = U_inf*sin(alpha);

Vn_inf = u_inf*nx + v_inf*ny;

% No penetration:
% A*lambda + Vn_inf = 0
lambda = A\(-Vn_inf(:));

%% 5. Source conservation
total_source = sum(lambda(:).*S(:));
source_residual = abs(total_source);

%% 6. Surface tangential velocity and Cp
Vt = zeros(N,1);

for i = 1:N

    uTotal = u_inf;
    vTotal = v_inf;

    for j = 1:N

        dx = xc(i)-x1(j);
        dy = yc(i)-y1(j);

        X =  dx*cos(beta(j)) + dy*sin(beta(j));
        Y = -dx*sin(beta(j)) + dy*cos(beta(j));

        if i == j
            uLocal = 0;
            vLocal = -0.5;
        else
            r1sq = X^2 + Y^2;
            r2sq = (X-S(j))^2 + Y^2;

            uLocal = (1/(4*pi))*log(r1sq/r2sq);

            theta1 = atan2(Y,X);
            theta2 = atan2(Y,X-S(j));

            vLocal = (1/(2*pi))*(theta2-theta1);
        end

        uTotal = uTotal + lambda(j)* ...
            (uLocal*cos(beta(j)) - vLocal*sin(beta(j)));

        vTotal = vTotal + lambda(j)* ...
            (uLocal*sin(beta(j)) + vLocal*cos(beta(j)));
    end

    Vt(i) = uTotal*tx(i) + vTotal*ty(i);
end

Cp_panel = 1 - (Vt/U_inf).^2;

theta_cp = atan2(yc,xc);
theta_cp(theta_cp < 0) = theta_cp(theta_cp < 0) + 2*pi;

%% 7. Analytical cylinder solution
theta_theory = linspace(0,2*pi,1000);
Cp_theory = 1 - 4*sin(theta_theory).^2;

Cp_exact_at_panels = 1 - 4*sin(theta_cp).^2;
Cp_error = Cp_panel - Cp_exact_at_panels;

Cp_min = min(Cp_panel);
Cp_max = max(Cp_panel);
RMSE_Cp = sqrt(mean(Cp_error.^2));
MAE_Cp = mean(abs(Cp_error));

%% 8. Print validation
fprintf('------------------------------------------------------------\n');
fprintf('SOURCE-PANEL VALIDATION\n');
fprintf('------------------------------------------------------------\n');
fprintf('Maximum source strength       = %+ .8f\n',max(lambda));
fprintf('Minimum source strength       = %+ .8f\n',min(lambda));
fprintf('Total source strength         = %+ .8e\n',total_source);
fprintf('Source conservation residual  = %.8e\n',source_residual);
fprintf('\n');
fprintf('Numerical Cp minimum          = %+ .8f\n',Cp_min);
fprintf('Numerical Cp maximum          = %+ .8f\n',Cp_max);
fprintf('Theoretical Cp minimum        = %+ .8f\n',-3);
fprintf('Theoretical Cp maximum        = %+ .8f\n',1);
fprintf('Cp RMSE                       = %.8e\n',RMSE_Cp);
fprintf('Cp MAE                        = %.8e\n',MAE_Cp);
fprintf('------------------------------------------------------------\n');

%% 9. Velocity field for streamlines
%
% IMPORTANT:
% The streamlines are calculated from the actual panel-induced velocity
% field. No analytical cylinder velocity is substituted here.

xg = linspace(-5,5,301);
yg = linspace(-5,5,301);
[Xg,Yg] = meshgrid(xg,yg);

Ug = zeros(size(Xg));
Vg = zeros(size(Yg));

for k = 1:numel(Xg)

    xp = Xg(k);
    yp = Yg(k);

    if hypot(xp,yp) <= R
        Ug(k) = NaN;
        Vg(k) = NaN;
        continue;
    end

    uTotal = u_inf;
    vTotal = v_inf;

    for j = 1:N

        dx = xp-x1(j);
        dy = yp-y1(j);

        X =  dx*cos(beta(j)) + dy*sin(beta(j));
        Y = -dx*sin(beta(j)) + dy*cos(beta(j));

        r1sq = X^2 + Y^2;
        r2sq = (X-S(j))^2 + Y^2;

        uLocal = (1/(4*pi))*log(r1sq/r2sq);

        theta1 = atan2(Y,X);
        theta2 = atan2(Y,X-S(j));

        vLocal = (1/(2*pi))*(theta2-theta1);

        uTotal = uTotal + lambda(j)* ...
            (uLocal*cos(beta(j)) - vLocal*sin(beta(j)));

        vTotal = vTotal + lambda(j)* ...
            (uLocal*sin(beta(j)) + vLocal*cos(beta(j)));
    end

    Ug(k) = uTotal;
    Vg(k) = vTotal;
end

%% 10. Task 4(a) - streamline figure
figure('Color','w','Position',[100 80 950 760]);
hold on;

% Use explicit seeds so the bending around the cylinder is clearly shown.
seedY = [-4.5 -4.0 -3.5 -3.0 -2.5 -2.0 -1.5 -1.20 -1.05 ...
         -0.90 -0.70 -0.50 -0.25 0 0.25 0.50 0.70 0.90 1.05 ...
          1.20 1.5 2.0 2.5 3.0 3.5 4.0 4.5];

seedX = -4.9*ones(size(seedY));

% Trace individual streamlines through the numerical velocity field.
for m = 1:numel(seedY)
    h = streamline(Xg,Yg,Ug,Vg,seedX(m),seedY(m));
    if ~isempty(h)
        set(h,'Color',[0 0.4470 0.7410],'LineWidth',1.0);
    end
end

% Cylinder surface
plot(R*cos(theta_nodes),R*sin(theta_nodes), ...
     'k-','LineWidth',2.2);

% Control points
plot(xc,yc,'k.','MarkerSize',5);

xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
title('Task 4(a): 64-Panel Source-Panel Flow over a Cylinder', ...
      'Interpreter','latex');

axis equal;
xlim([-5 5]);
ylim([-5 5]);
grid on;
box on;
set(gca,'FontSize',11);

% Create clean legend handles
hStream = plot(nan,nan,'-','Color',[0 0.4470 0.7410],'LineWidth',1.0);
hCyl = plot(nan,nan,'k-','LineWidth',2.2);
hCP = plot(nan,nan,'k.','MarkerSize',12);

legend([hStream hCyl hCP], ...
       {'Streamlines','Cylinder','Control points'}, ...
       'Location','northeast');

%% 11. Task 4(b) - Cp comparison
figure('Color','w','Position',[120 100 1000 650]);
hold on;

plot(theta_cp,Cp_panel, ...
     'ko','MarkerSize',4.5,'MarkerFaceColor','k');

plot(theta_theory,Cp_theory, ...
     'b-','LineWidth',1.8);

xlabel('$\theta$ [rad]','Interpreter','latex');
ylabel('$C_p$','Interpreter','latex');

title('Task 4(b): Cylinder Pressure Coefficient Comparison', ...
      'Interpreter','latex');

legend({'64-panel source method','$C_p=1-4\sin^2\theta$'}, ...
       'Interpreter','latex','Location','best');

xlim([0 2*pi]);
ylim([-3.1 1.1]);

xticks([0 pi/2 pi 3*pi/2 2*pi]);
xticklabels({'0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'});

grid on;
box on;
set(gca,'FontSize',11);

%% 12. Final output
fprintf('\n============================================================\n');
fprintf(' FINAL TASK 4(a)/(b) RESULTS\n');
fprintf('============================================================\n');
fprintf('Panels                    : %d\n',N);
fprintf('Cylinder radius           : %.4f\n',R);
fprintf('Freestream velocity       : %.4f\n',U_inf);
fprintf('Cp minimum                : %.6f\n',Cp_min);
fprintf('Cp maximum                : %.6f\n',Cp_max);
fprintf('Cp RMSE                   : %.6e\n',RMSE_Cp);
fprintf('Cp MAE                    : %.6e\n',MAE_Cp);
fprintf('Source residual           : %.6e\n',source_residual);
fprintf('============================================================\n');
fprintf(' Task 4(a) and 4(b) COMPLETE.\n');
fprintf('============================================================\n');
