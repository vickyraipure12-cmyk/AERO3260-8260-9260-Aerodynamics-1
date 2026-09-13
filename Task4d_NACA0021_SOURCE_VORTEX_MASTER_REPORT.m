%% AERO3260 Assignment 1 - Task 4(d) MASTER / REPORT VERSION
% NACA 0021, 64 panels, Uinf = 1, alpha = 10 deg
% Constant-strength source panels + one constant-strength vortex sheet
% with Kutta condition at the trailing edge.
%
% This version fixes:
%   1) row/column dimension errors
%   2) Kutta residual calculation
%   3) Cp plotting (upper/lower surfaces plotted separately)
%   4) force-coefficient sign convention
%   5) x-axis label
%   6) cleaner streamline and trailing-edge streamline plots
%   7) numerical data export
%
% NOTE:
% The NACA geometry is closed at the trailing edge (same convention as
% the corrected Task 4(c) MATLAB source-panel calculation).

clear; clc; close all;

%% ================================================================
% 1. ASSIGNMENT INPUTS
% ================================================================
N = 64;                    % number of panels
Uinf = 1.0;                % freestream velocity
alpha_deg = 10.0;          % angle of attack [deg]
alpha = deg2rad(alpha_deg);
t = 0.21;                  % NACA 0021 thickness ratio
c = 1.0;                   % chord

%% ================================================================
% 2. NACA 0021 GEOMETRY
% ================================================================
nHalf = N/2 + 1;
beta = linspace(0,pi,nHalf);
xHalf = (1-cos(beta))/2;

yt = 5*t*( ...
      0.2969*sqrt(xHalf) ...
    - 0.1260*xHalf ...
    - 0.3516*xHalf.^2 ...
    + 0.2843*xHalf.^3 ...
    - 0.1015*xHalf.^4);

% Close trailing edge to match Task 4(c) MATLAB geometry
yt(end) = 0;

% Boundary ordering:
% upper trailing edge -> leading edge -> lower trailing edge
x = [fliplr(xHalf), xHalf(2:end)];
y = [fliplr(yt),    -yt(2:end)];

% Force all geometry arrays to columns
x = x(:);
y = y(:);

assert(numel(x)==N+1,'Geometry must contain N+1 nodes.');

%% ================================================================
% 3. PANEL GEOMETRY
% ================================================================
x1 = x(1:N);
y1 = y(1:N);
x2 = x(2:N+1);
y2 = y(2:N+1);

dx = x2-x1;
dy = y2-y1;

S = hypot(dx,dy);                    % panel lengths
phi = atan2(dy,dx);                  % panel angles

xc = (x1+x2)/2;                      % control points
yc = (y1+y2)/2;

% For clockwise boundary ordering:
% n_in = (-sin(phi), cos(phi)) points INTO the airfoil.
nx = -sin(phi);
ny =  cos(phi);

% Panel tangent follows the panel ordering.
tx = cos(phi);
ty = sin(phi);

%% ================================================================
% 4. SOURCE/VORTEX INFLUENCE MATRICES
% ================================================================
As = zeros(N,N);       % source -> normal velocity
Av = zeros(N,1);       % common vortex -> normal velocity

Ts = zeros(N,N);       % source -> tangential velocity
Tv = zeros(N,1);       % common vortex -> tangential velocity

for i = 1:N
    for j = 1:N

        % Control point i in local coordinates of panel j
        X = (xc(i)-x1(j))*cos(phi(j)) ...
          + (yc(i)-y1(j))*sin(phi(j));

        Y = -(xc(i)-x1(j))*sin(phi(j)) ...
          + (yc(i)-y1(j))*cos(phi(j));

        if i == j

            % Correct self influence for the source panel
            uLocal = 0.0;
            vLocal = -0.5;

        else

            r1sq = X^2 + Y^2;
            r2sq = (X-S(j))^2 + Y^2;

            uLocal = (1/(4*pi))*log(r1sq/r2sq);

            dTheta = atan2(Y,X-S(j)) - atan2(Y,X);

            % Wrap angle difference to [-pi,pi] to avoid branch jumps.
            dTheta = atan2(sin(dTheta),cos(dTheta));

            vLocal = (1/(2*pi))*dTheta;
        end

        % ---------------- SOURCE PANEL ----------------
        uS = uLocal*cos(phi(j)) - vLocal*sin(phi(j));
        vS = uLocal*sin(phi(j)) + vLocal*cos(phi(j));

        As(i,j) = uS*nx(i) + vS*ny(i);
        Ts(i,j) = uS*tx(i) + vS*ty(i);

        % ---------------- VORTEX PANEL ----------------
        % A unit vortex panel is a 90-degree rotation of the
        % corresponding source-panel velocity field.
        uVL = -vLocal;
        vVL =  uLocal;

        uV = uVL*cos(phi(j)) - vVL*sin(phi(j));
        vV = uVL*sin(phi(j)) + vVL*cos(phi(j));

        Av(i) = Av(i) + uV*nx(i) + vV*ny(i);
        Tv(i) = Tv(i) + uV*tx(i) + vV*ty(i);

    end
end

%% ================================================================
% 5. FREESTREAM COMPONENTS
% ================================================================
uInf = Uinf*cos(alpha);
vInf = Uinf*sin(alpha);

VnInf = uInf*nx + vInf*ny;
VtInf = uInf*tx + vInf*ty;

%% ================================================================
% 6. 65 x 65 HESS-SMITH SYSTEM
% ================================================================
% Unknown vector:
% [lambda_1 ... lambda_64 gamma]^T
%
% Equations:
% 64 no-penetration equations
% 1 Kutta equation

A = zeros(N+1,N+1);
b = zeros(N+1,1);

% No penetration
A(1:N,1:N) = As;
A(1:N,N+1) = Av;
b(1:N) = -VnInf;

% Kutta condition.
% The first panel is immediately below/above the trailing edge according
% to the clockwise ordering, while panel N is the other TE panel.
% Their local tangents point in opposite directions, so smooth flow
% requires Vt(1) + Vt(N) = 0.
A(N+1,1:N) = Ts(1,:) + Ts(N,:);
A(N+1,N+1) = Tv(1) + Tv(N);
b(N+1) = -(VtInf(1) + VtInf(N));

%% ================================================================
% 7. SOLVE
% ================================================================
q = A\b;

lambda = q(1:N);
gamma = q(N+1);

%% ================================================================
% 8. SURFACE VELOCITY AND PRESSURE COEFFICIENT
% ================================================================
Vt = VtInf + Ts*lambda + Tv*gamma;

Cp = 1 - (Vt/Uinf).^2;

assert(iscolumn(Cp),'Cp must be a column vector.');
assert(numel(Cp)==N,'Cp must contain one value per panel.');

%% ================================================================
% 9. AERODYNAMIC COEFFICIENTS
% ================================================================
% nx,ny defined above point INTO the airfoil for this clockwise boundary.
% Therefore the physical outward normal is:
% n_out = -[nx,ny]
%
% Pressure force:
% dF = -Cp*n_out*dS = Cp*[nx,ny]*dS
%
% Body-axis force coefficients:
Cx = sum(Cp.*nx.*S/c);
Cy = sum(Cp.*ny.*S/c);
%
% Transform body-axis x-y forces to lift/drag relative to freestream:
Cd =  Cx*cos(alpha) + Cy*sin(alpha);
Cl = -Cx*sin(alpha) + Cy*cos(alpha);

% Kutta and source-conservation checks
kuttaResidual = Vt(1) + Vt(N);
sourceResidual = sum(lambda.*S);

% Total circulation of the constant vortex sheet
Gamma = gamma*sum(S);

%% ================================================================
% 10. PRINT RESULTS
% ================================================================
fprintf('\n============================================================\n');
fprintf('AERO3260 TASK 4(d) - NACA 0021 SOURCE + VORTEX PANEL METHOD\n');
fprintf('============================================================\n');
fprintf('Number of panels          = %d\n',N);
fprintf('Angle of attack           = %.2f deg\n',alpha_deg);
fprintf('Freestream velocity       = %.4f\n',Uinf);
fprintf('Vortex strength gamma     = %.8f\n',gamma);
fprintf('Total circulation Gamma   = %.8f\n',Gamma);
fprintf('Source conservation       = %.6e\n',sourceResidual);
fprintf('Kutta residual             = %.6e\n',kuttaResidual);
fprintf('CL                         = %.8f\n',Cl);
fprintf('CD                         = %.8f\n',Cd);
fprintf('Cx                         = %.8f\n',Cx);
fprintf('Cy                         = %.8f\n',Cy);
fprintf('Cp minimum                 = %.8f\n',min(Cp));
fprintf('Cp maximum                 = %.8f\n',max(Cp));
fprintf('Condition number           = %.4e\n',cond(A));
fprintf('============================================================\n');

%% ================================================================
% 11. FIGURE 1 - AIRFOIL GEOMETRY
% ================================================================
figure('Color','w','Name','Task 4(d) - Geometry');

plot(x,y,'k-','LineWidth',1.8);
hold on;
plot(xc,yc,'ko','MarkerSize',3.5,'MarkerFaceColor','k');

grid on;
axis equal;
xlim([-0.05 1.05]);
ylim([-0.25 0.25]);

xlabel('x/c','FontSize',12);
ylabel('y/c','FontSize',12);
title('NACA 0021 Geometry and Panel Control Points, \alpha = 10^\circ',...
    'FontSize',14);

set(gca,'FontSize',11,'LineWidth',1);

%% ================================================================
% 12. FIGURE 2 - CLEAN Cp DISTRIBUTION
% ================================================================
% First 32 panels = upper surface (TE -> LE)
% Last 32 panels = lower surface (LE -> TE)

nSide = N/2;

xu = xc(1:nSide);
Cpu = Cp(1:nSide);

xl = xc(nSide+1:end);
Cpl = Cp(nSide+1:end);

% Sort each surface from LE -> TE for a clean x/c plot
[xu,iu] = sort(xu);
Cpu = Cpu(iu);

[xl,il] = sort(xl);
Cpl = Cpl(il);

figure('Color','w','Name','Task 4(d) - Cp');

plot(xu,Cpu,'k-o','LineWidth',1.4,'MarkerSize',4,...
    'MarkerFaceColor','k');
hold on;
plot(xl,Cpl,'k-s','LineWidth',1.4,'MarkerSize',4,...
    'MarkerFaceColor','w');

set(gca,'YDir','reverse');
grid on;
box on;

xlabel('x/c','FontSize',12);
ylabel('C_p','FontSize',12);
title('NACA 0021 Surface Pressure Distribution, \alpha = 10^\circ',...
    'FontSize',14);

legend('Upper surface','Lower surface','Location','best');
xlim([0 1]);

set(gca,'FontSize',11,'LineWidth',1);

%% ================================================================
% 13. STREAMLINE VELOCITY FIELD
% ================================================================
xMin = -1.5;
xMax =  2.5;
yMin = -1.5;
yMax =  1.5;

ng = 220;

[Xg,Yg] = meshgrid(linspace(xMin,xMax,ng),...
                   linspace(yMin,yMax,ng));

[Ug,Vg] = panelVelocity(Xg,Yg,...
    x1,y1,S,phi,lambda,gamma,Uinf,alpha);

% Remove interior points
inside = inpolygon(Xg,Yg,x,y);
Ug(inside) = NaN;
Vg(inside) = NaN;

%% ================================================================
% 14. FIGURE 3 - STREAMLINES + TRAILING-EDGE STREAMLINE
% ================================================================
figure('Color','w','Name','Task 4(d) - Streamlines');

% General flow streamlines
streamslice(Xg,Yg,Ug,Vg,1.6);
hold on;

% Airfoil body
fill(x,y,'w','EdgeColor','k','LineWidth',1.8);

% Trailing-edge streamline.
% Do NOT start exactly at (1,0), because the ideal sharp TE is a
% mathematical singular point. Start immediately downstream.
xSeed = 1.002;
ySeed = 0.0005;

hTE = streamline(Xg,Yg,Ug,Vg,xSeed,ySeed);

if ~isempty(hTE)
    set(hTE,'LineWidth',2.5);
end

% Mark trailing edge
plot(1,0,'ko','MarkerFaceColor','k','MarkerSize',6);

grid on;
box on;
axis equal;

xlabel('x/c','FontSize',12);
ylabel('y/c','FontSize',12);
title('NACA 0021 Streamlines, \alpha = 10^\circ',...
    'FontSize',14);

xlim([xMin xMax]);
ylim([yMin yMax]);

set(gca,'FontSize',11,'LineWidth',1);

%% ================================================================
% 15. SAVE PANEL DATA
% ================================================================
% Columns:
% 1 x/c
% 2 y/c
% 3 Cp
% 4 Vt/Uinf
% 5 source strength lambda
% 6 panel length
panelResults = [xc/c, yc/c, Cp, Vt/Uinf, lambda, S];

writematrix(panelResults,...
    'Task4d_NACA0021_source_vortex_panel_results.csv');

fprintf('\nSaved numerical results:\n');
fprintf('Task4d_NACA0021_source_vortex_panel_results.csv\n\n');

%% ================================================================
% LOCAL FUNCTION - VELOCITY FIELD
% ================================================================
function [u,v] = panelVelocity(Xp,Yp,...
    x1,y1,S,phi,lambda,gamma,Uinf,alpha)

u = Uinf*cos(alpha)*ones(size(Xp));
v = Uinf*sin(alpha)*ones(size(Xp));

for j = 1:length(S)

    X = (Xp-x1(j))*cos(phi(j)) ...
      + (Yp-y1(j))*sin(phi(j));

    Y = -(Xp-x1(j))*sin(phi(j)) ...
      + (Yp-y1(j))*cos(phi(j));

    r1sq = max(X.^2 + Y.^2,1e-12);
    r2sq = max((X-S(j)).^2 + Y.^2,1e-12);

    uLocal = (1/(4*pi))*log(r1sq./r2sq);

    dTheta = atan2(Y,X-S(j))-atan2(Y,X);
    dTheta = atan2(sin(dTheta),cos(dTheta));

    vLocal = (1/(2*pi))*dTheta;

    % Source panel velocity
    uS = uLocal*cos(phi(j))-vLocal*sin(phi(j));
    vS = uLocal*sin(phi(j))+vLocal*cos(phi(j));

    % Vortex panel velocity
    uVL = -vLocal;
    vVL =  uLocal;

    uV = uVL*cos(phi(j))-vVL*sin(phi(j));
    vV = uVL*sin(phi(j))+vVL*cos(phi(j));

    u = u + lambda(j)*uS + gamma*uV;
    v = v + lambda(j)*vS + gamma*vV;
end
end
