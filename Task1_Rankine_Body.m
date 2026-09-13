%% AERO3260 Assignment 1 - Task 1 (FINAL VERIFIED)
% Rankine body formed by a uniform stream + source/sink pair.
% Student ID: 550361509
%
% Task 1(b): Uinf = sigma = x0 = 1
% Task 1(c): L = 50 m, T = 10 m, Uinf = 5 m/s
%
% IMPORTANT:
% atan2 is used so that the source/sink streamfunction has the correct
% quadrant/branch behaviour.

clear; close all; clc

%% ================================================================
% TASK 1(b) - Normalised Rankine body
% ================================================================
Uinf = 1;
sigma = 1;
x0 = 1;

x = linspace(-3,3,1200);
y = linspace(-2.5,2.5,1000);
[X,Y] = meshgrid(x,y);

psi = Uinf.*Y + sigma/(2*pi) .* ...
    (atan2(Y,X+x0) - atan2(Y,X-x0));

% Stagnation points are obtained from u = 0 on y = 0:
% xs^2 = x0^2 + sigma*x0/(pi*Uinf)
xs = sqrt(x0^2 + sigma*x0/(pi*Uinf));

fprintf('\n================ TASK 1(b) ================\n');
fprintf('Uinf = %.4f, sigma = %.4f, x0 = %.4f\n',Uinf,sigma,x0);
fprintf('Stagnation points: x = +/- %.8f\n',xs);

figure('Color','w');
contour(X,Y,psi,45,'LineWidth',0.8); hold on
plot([-xs xs],[0 0],'ko','MarkerFaceColor','k','MarkerSize',6)

% Highlight the Rankine-body boundary (psi = 0)
contour(X,Y,psi,[0 0],'k','LineWidth',2.2);

axis equal
xlim([-3 3]); ylim([-2.5 2.5]);
xlabel('x','FontSize',12)
ylabel('y','FontSize',12)
title('Rankine Body Streamlines: U_\infty = \sigma = x_0 = 1')
grid on; box on

%% ================================================================
% TASK 1(c) - Submarine
% ================================================================
L = 50;                 % m
T = 10;                 % m
Uinf = 5;               % m/s

% The two stagnation points are the bow and stern:
% 2*xs = L
xs_target = L/2;

% Two equations must be satisfied:
%
% (1) xs^2 = x0^2 + sigma*x0/(pi*Uinf)
%
% (2) The body streamline passes through (0,T/2):
%
% Uinf*T/2 + sigma/(2*pi) * ...
% [2*atan(T/(2*x0)) - pi] = 0
%
% Solve the two nonlinear equations simultaneously.

fun = @(q) [ ...
    q(1)^2 + q(2)*q(1)/(pi*Uinf) - xs_target^2;
    Uinf*T/2 + q(2)/(2*pi) * ...
    (2*atan(T/(2*q(1))) - pi) ];

q0 = [23.2 57.8];
opts = optimoptions('fsolve','Display','off',...
    'FunctionTolerance',1e-12,'StepTolerance',1e-12);

q = fsolve(fun,q0,opts);
x0_sub = q(1);
sigma_sub = q(2);

% Verify residuals independently.
residuals = fun([x0_sub sigma_sub]);
xs_check = sqrt(x0_sub^2 + sigma_sub*x0_sub/(pi*Uinf));

fprintf('\n================ TASK 1(c) ================\n');
fprintf('Required submarine length  L = %.4f m\n',L);
fprintf('Required thickness         T = %.4f m\n',T);
fprintf('Cruise speed          Uinf = %.4f m/s\n',Uinf);
fprintf('Stagnation location    xs = %.8f m\n',xs_target);
fprintf('Source/sink x0             = %.8f m\n',x0_sub);
fprintf('Source/sink strength sigma = %.8f m^2/s\n',sigma_sub);
fprintf('Source/sink spacing 2*x0   = %.8f m\n',2*x0_sub);
fprintf('Stagnation check            = %.8f m\n',xs_check);
fprintf('Equation residual 1         = %.3e\n',residuals(1));
fprintf('Equation residual 2         = %.3e\n',residuals(2));

%% Reproduce the Rankine-body plot using the calculated submarine values
x = linspace(-26,26,1400);
y = linspace(-8,8,1000);
[X,Y] = meshgrid(x,y);

psi_sub = Uinf.*Y + sigma_sub/(2*pi) .* ...
    (atan2(Y,X+x0_sub) - atan2(Y,X-x0_sub));

figure('Color','w');
contour(X,Y,psi_sub,45,'LineWidth',0.8); hold on
contour(X,Y,psi_sub,[0 0],'k','LineWidth',2.5);
plot([-xs_target xs_target],[0 0],'ko',...
    'MarkerFaceColor','k','MarkerSize',6)

% Mark the specified maximum-thickness points.
plot(0,T/2,'ks','MarkerFaceColor','w','MarkerSize',7)
plot(0,-T/2,'ks','MarkerFaceColor','w','MarkerSize',7)

axis equal
xlim([-26 26]); ylim([-8 8]);
xlabel('x (m)','FontSize',12)
ylabel('y (m)','FontSize',12)
title('Rankine-Body Representation of the 50 m Submarine')
legend('Streamlines','Rankine body','Stagnation points',...
       'Thickness points','Location','best')
grid on; box on

% Final direct thickness check:
psi_mid = Uinf*(T/2) + sigma_sub/(2*pi) * ...
    (atan2(T/2,x0_sub)-atan2(T/2,-x0_sub));

fprintf('Direct body-streamline check at (0,T/2): psi = %.3e\n',psi_mid);
fprintf('============================================\n\n');
