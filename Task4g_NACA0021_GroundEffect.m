%% ========================================================================
% AERO8260 / AERO9260 - Aerodynamics 1
% Assignment 1 - Task 4(g)
%
% NACA 0021 Aerofoil in Ground Effect
% alpha = 10 deg
% Leading edge height = 0.5 chord
%
% Method:
%   1. 64-panel constant-source + constant-vortex panel method
%   2. Kutta condition at the trailing edge
%   3. Ground-plane image system
%      - source image: same strength
%      - vortex image: opposite strength
%   4. Cp field and streamlines are evaluated over the flow domain
%
% Output:
%   Task4g_NACA0021_GroundEffect.png
%
% Note:
%   This is an inviscid potential-flow image-system ground-effect
%   visualization. Viscous boundary-layer, transition, separation and
%   wake-loss effects are not included.
% ========================================================================

clear;
clc;
close all;

%% ------------------------ INPUT PARAMETERS -----------------------------

N = 64;                         % Number of panels
alpha_deg = 10;                 % Angle of attack [deg]
alpha = deg2rad(alpha_deg);

h_LE = 0.50;                    % Leading-edge height / chord
Uinf = 1.0;                     % Freestream velocity
c = 1.0;                        % Chord

t = 0.21;                       % NACA 0021 thickness ratio

%% ------------------------ NACA 0021 GEOMETRY ----------------------------

% Cosine-spaced x coordinates
beta = linspace(0, pi, N/2 + 1);
x = (1 - cos(beta))/2;

% NACA 00xx thickness distribution
yt = 5*t*( ...
      0.2969*sqrt(x) ...
    - 0.1260*x ...
    - 0.3516*x.^2 ...
    + 0.2843*x.^3 ...
    - 0.1015*x.^4 );

% Upper surface: trailing edge -> leading edge
xu = fliplr(x);
yu = fliplr(yt);

% Lower surface: leading edge -> trailing edge
xl = x(2:end);
yl = -yt(2:end);

% Closed aerofoil
x_airfoil = [xu, xl];
y_airfoil = [yu, yl];

% Place leading edge 0.5 chords above ground
y_airfoil = y_airfoil + h_LE;

N_actual = length(x_airfoil) - 1;

%% ------------------------ PANEL GEOMETRY -------------------------------

XA = x_airfoil(1:end-1);
YA = y_airfoil(1:end-1);

XB = x_airfoil(2:end);
YB = y_airfoil(2:end);

dx = XB - XA;
dy = YB - YA;

S = sqrt(dx.^2 + dy.^2);

% Panel tangent vectors
tx = dx ./ S;
ty = dy ./ S;

% Outward normals for clockwise boundary
nx = ty;
ny = -tx;

% Control points
XC = 0.5*(XA + XB);
YC = 0.5*(YA + YB);

%% ------------------------ FREESTREAM -----------------------------------

Vx_inf = Uinf*cos(alpha);
Vy_inf = Uinf*sin(alpha);

%% ------------------------ INFLUENCE MATRICES ---------------------------

As = zeros(N_actual, N_actual);      % Normal source influence
Av = zeros(N_actual, N_actual);      % Normal vortex influence
Bs = zeros(N_actual, N_actual);      % Tangential source influence
Bv = zeros(N_actual, N_actual);      % Tangential vortex influence

for i = 1:N_actual

    for j = 1:N_actual

        % Coordinates of control point i relative to panel j
        rx = XC(i) - XA(j);
        ry = YC(i) - YA(j);

        % Transform into panel-local coordinates
        xloc = rx*tx(j) + ry*ty(j);
        yloc = rx*nx(j) + ry*ny(j);

        Sj = S(j);

        if i == j

            % Principal-value/self influence
            us = 0.0;
            vs = 0.5;

            uv = 0.0;
            vv = 0.0;

        else

            r1sq = xloc^2 + yloc^2;
            r2sq = (xloc-Sj)^2 + yloc^2;

            r1sq = max(r1sq, 1e-14);
            r2sq = max(r2sq, 1e-14);

            % Constant-strength source panel
            us = (1/(4*pi))*log(r1sq/r2sq);

            theta1 = atan2(yloc, xloc);
            theta2 = atan2(yloc, xloc-Sj);

            vs = (theta2-theta1)/(2*pi);

            % Constant-strength vortex panel
            uv = -vs;
            vv = us;

        end

        % Convert source velocity to global coordinates
        u_source = us*tx(j) + vs*nx(j);
        v_source = us*ty(j) + vs*ny(j);

        % Convert vortex velocity to global coordinates
        u_vortex = uv*tx(j) + vv*nx(j);
        v_vortex = uv*ty(j) + vv*ny(j);

        % Normal influence
        As(i,j) = u_source*nx(i) + v_source*ny(i);
        Av(i,j) = u_vortex*nx(i) + v_vortex*ny(i);

        % Tangential influence
        Bs(i,j) = u_source*tx(i) + v_source*ty(i);
        Bv(i,j) = u_vortex*tx(i) + v_vortex*ty(i);

    end
end

%% ------------------------ SOLVE PANEL SYSTEM ---------------------------

% No-penetration boundary condition
RHS = -(Vx_inf*nx + Vy_inf*ny)';

% Unknowns:
%   sigma(1:N_actual) = source strengths
%   Gamma              = constant vortex-sheet strength
A_system = zeros(N_actual+1, N_actual+1);
b_system = zeros(N_actual+1, 1);

A_system(1:N_actual,1:N_actual) = As;
A_system(1:N_actual,N_actual+1) = sum(Av,2);

b_system(1:N_actual) = RHS;

% Kutta condition:
% Vt(upper TE) + Vt(lower TE) = 0
A_system(N_actual+1,1:N_actual) = Bs(1,:) + Bs(end,:);
A_system(N_actual+1,N_actual+1) = ...
    sum(Bv(1,:)) + sum(Bv(end,:));

b_system(N_actual+1) = ...
    -(Vx_inf*tx(1) + Vy_inf*ty(1)) ...
    -(Vx_inf*tx(end) + Vy_inf*ty(end));

% Solve for source strengths and vortex-sheet strength
solution = A_system \ b_system;

sigma = solution(1:N_actual);
Gamma = solution(N_actual+1);

%% ------------------------ SURFACE VELOCITY / Cp ------------------------

Vt = zeros(N_actual,1);

for i = 1:N_actual

    Vt(i) = ...
        Vx_inf*tx(i) + Vy_inf*ty(i) ...
        + Bs(i,:)*sigma ...
        + sum(Bv(i,:))*Gamma;

end

Cp_surface = 1 - (Vt/Uinf).^2;

%% ------------------------ FLOW-FIELD GRID ------------------------------

x_min = -1.5;
x_max =  2.8;
y_min =  0.0;
y_max =  1.95;

nx_field = 220;
ny_field = 150;

xg = linspace(x_min, x_max, nx_field);
yg = linspace(y_min, y_max, ny_field);

[X,Y] = meshgrid(xg,yg);

% Freestream
U = Vx_inf*ones(size(X));
V = Vy_inf*ones(size(X));

%% ------------------------ PHYSICAL AEROFOIL ----------------------------

for j = 1:N_actual

    [u_panel,v_panel] = panelVelocityField( ...
        X,Y, ...
        XA(j),YA(j), ...
        XB(j),YB(j), ...
        sigma(j),Gamma);

    U = U + u_panel;
    V = V + v_panel;

end

%% ------------------------ GROUND IMAGE SYSTEM --------------------------

% Reflect the aerofoil about y = 0.
%
% Source image:
%       same sign
%
% Vortex image:
%       opposite sign
%
% This provides the ground-plane image-system representation.

for j = 1:N_actual

    XA_img = XA(j);
    YA_img = -YA(j);

    XB_img = XB(j);
    YB_img = -YB(j);

    % Source image
    [u_source_img,v_source_img] = panelVelocityField( ...
        X,Y, ...
        XA_img,YA_img, ...
        XB_img,YB_img, ...
        sigma(j),0);

    % Vortex image
    [u_vortex_img,v_vortex_img] = panelVelocityField( ...
        X,Y, ...
        XA_img,YA_img, ...
        XB_img,YB_img, ...
        0,-Gamma);

    U = U + u_source_img + u_vortex_img;
    V = V + v_source_img + v_vortex_img;

end

%% ------------------------ Cp FIELD -------------------------------------

Vmag = sqrt(U.^2 + V.^2);

Cp = 1 - (Vmag/Uinf).^2;

% Mask aerofoil interior
inside = inpolygon(X,Y,x_airfoil,y_airfoil);

Cp(inside) = NaN;
U(inside) = NaN;
V(inside) = NaN;

%% ------------------------ FIGURE ---------------------------------------

figure( ...
    'Color','w', ...
    'Position',[100 80 1250 680]);

% Pressure coefficient field
contourf(X,Y,Cp,30,'LineColor','none');
hold on;

% Fixed colour limits for consistent presentation
caxis([-5 1]);

%% ------------------------ STREAMLINES ----------------------------------

% Use a limited number of starting points so that the Cp field remains
% clearly visible. Starting from the upstream boundary gives clean,
% physically interpretable streamlines.

n_stream = 22;

start_x = x_min*ones(1,n_stream);
start_y = linspace(0.08,1.86,n_stream);

for k = 1:n_stream

    x0 = start_x(k);
    y0 = start_y(k);

    if inpolygon(x0,y0,x_airfoil,y_airfoil)
        continue;
    end

    try

        stream = stream2(X,Y,U,V,x0,y0);

        if ~isempty(stream) && ~isempty(stream{1})

            plot( ...
                stream{1}(:,1), ...
                stream{1}(:,2), ...
                'k-', ...
                'LineWidth',0.55);

        end

    catch
        % Skip unsuccessful streamline integrations
    end

end

%% ------------------------ TRAILING-EDGE STREAMLINE ---------------------

% The NACA geometry is defined TE -> LE on the upper surface and
% LE -> TE on the lower surface. Therefore the first point is the TE.

x_TE = x_airfoil(1);
y_TE = y_airfoil(1);

% Start just downstream of the TE to avoid evaluating exactly at the
% panel singularity.
x_start_TE = x_TE + 0.015;
y_start_TE = y_TE;

try

    TE_stream = stream2( ...
        X,Y,U,V, ...
        x_start_TE,y_start_TE);

    if ~isempty(TE_stream) && ~isempty(TE_stream{1})

        TE_xy = TE_stream{1};

        % Highlight the trailing-edge streamline
        plot( ...
            TE_xy(:,1), ...
            TE_xy(:,2), ...
            '-', ...
            'Color',[1 0.85 0], ...
            'LineWidth',2.5);

        % Label it without covering the main flow field
        label_id = min(35,size(TE_xy,1));

        text( ...
            TE_xy(label_id,1), ...
            TE_xy(label_id,2)+0.055, ...
            'Trailing-edge streamline', ...
            'FontSize',9, ...
            'Color','k', ...
            'BackgroundColor','w', ...
            'Margin',2, ...
            'EdgeColor',[0.5 0.5 0.5]);

    end

catch
    % Continue if TE streamline cannot be integrated
end

%% ------------------------ AEROFOIL OUTLINE -----------------------------

plot( ...
    x_airfoil, ...
    y_airfoil, ...
    'k-', ...
    'LineWidth',2.0);

%% ------------------------ GROUND PLANE --------------------------------

plot( ...
    [x_min x_max], ...
    [0 0], ...
    'k-', ...
    'LineWidth',2.0);

% Ground-plane label
text( ...
    2.15,0.035, ...
    'Ground plane', ...
    'FontSize',9, ...
    'FontWeight','bold', ...
    'BackgroundColor','w', ...
    'Margin',1.5);

%% ------------------------ FORMATTING -----------------------------------

axis equal;
xlim([x_min x_max]);
ylim([y_min y_max]);

xlabel('x/c','FontSize',11);
ylabel('y/c','FontSize',11);

title( ...
    sprintf( ...
    'NACA 0021 in Ground Effect: Source--Vortex Panel Image System, \\alpha = %g^\\circ, h_{LE}/c = %.1f', ...
    alpha_deg,h_LE), ...
    'FontSize',12, ...
    'FontWeight','bold');

cb = colorbar;
cb.Label.String = 'C_p';
cb.FontSize = 10;

box on;
set(gca,'FontSize',10);

% Keep the figure clean and publication-like
set(gca,'Layer','top');

%% ------------------------ SAVE FIGURE ----------------------------------

% Export a high-resolution PNG for the report
exportgraphics( ...
    gcf, ...
    'Task4g_NACA0021_GroundEffect.png', ...
    'Resolution',300);

fprintf('\n');
fprintf('============================================================\n');
fprintf(' Task 4(g) completed successfully.\n');
fprintf(' NACA 0021, alpha = %.1f deg, h_LE/c = %.2f\n', ...
    alpha_deg,h_LE);
fprintf(' Number of panels = %d\n',N_actual);
fprintf(' Vortex-sheet strength Gamma = %.8f\n',Gamma);
fprintf(' Figure saved as:\n');
fprintf(' Task4g_NACA0021_GroundEffect.png\n');
fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTION: VELOCITY INDUCED BY ONE CONSTANT SOURCE/VORTEX PANEL
% ========================================================================

function [u,v] = panelVelocityField( ...
    X,Y,XA,YA,XB,YB,sigma,Gamma)

    % Panel geometry
    dx = XB-XA;
    dy = YB-YA;

    S = sqrt(dx^2 + dy^2);

    tx = dx/S;
    ty = dy/S;

    % Normal direction
    nx = ty;
    ny = -tx;

    % Coordinates relative to panel start point
    rx = X-XA;
    ry = Y-YA;

    % Transform into local panel coordinates
    xloc = rx*tx + ry*ty;
    yloc = rx*nx + ry*ny;

    % Squared distances to panel end points
    r1sq = xloc.^2 + yloc.^2;
    r2sq = (xloc-S).^2 + yloc.^2;

    % Protect against exact singularities
    r1sq = max(r1sq,1e-12);
    r2sq = max(r2sq,1e-12);

    % Source panel influence
    us = (1/(4*pi))*log(r1sq./r2sq);

    theta1 = atan2(yloc,xloc);
    theta2 = atan2(yloc,xloc-S);

    vs = (theta2-theta1)/(2*pi);

    % Vortex panel influence
    uv = -vs;
    vv = us;

    % Total local velocity
    u_local = sigma*us + Gamma*uv;
    v_local = sigma*vs + Gamma*vv;

    % Transform local velocity to global coordinates
    u = u_local*tx + v_local*nx;
    v = u_local*ty + v_local*ny;

    % Remove numerical non-finite values
    u(~isfinite(u)) = 0;
    v(~isfinite(v)) = 0;

end
