%% AERO3260 Assignment 1 - Task 4(c)
% NACA 0021 constant-source panel method at alpha = 0 deg
% 64 panels + comparison with inviscid XFOIL
%
% Required result:
%   - Symmetric Cp distribution at alpha = 0 deg
%   - Comparison with inviscid XFOIL
%
% XFOIL file expected in the same folder:
%   xfoil_NACA0021_alpha0_inviscid.txt

clear; clc; close all;

%% -------------------- INPUTS --------------------
N = 64;                 % number of panels
t = 0.21;               % NACA 0021 thickness ratio
alpha_deg = 0.0;        % angle of attack
Uinf = 1.0;              % freestream speed
rho = 1.0;              % density (only needed for force checks)

xfoilFile = 'xfoil_NACA0021_alpha0_inviscid.txt';

%% ---------------- NACA 0021 GEOMETRY ----------------
% Cosine spacing gives good resolution near LE and TE.
beta = linspace(0,pi,N/2+1);
x = (1-cos(beta))/2;

% NACA 00xx thickness distribution.
yt = 5*t*( ...
      0.2969*sqrt(x) ...
    - 0.1260*x ...
    - 0.3516*x.^2 ...
    + 0.2843*x.^3 ...
    - 0.1015*x.^4 );

% Build a closed contour:
% upper surface: TE -> LE
% lower surface: LE -> TE
% This ordering is counter-clockwise.
xUpper = fliplr(x);
yUpper = fliplr(yt);

xLower = x(2:end);
yLower = -yt(2:end);

X = [xUpper, xLower];
Y = [yUpper, yLower];

% Remove duplicate final point if present and close explicitly.
if abs(X(end)-X(1)) < 1e-14 && abs(Y(end)-Y(1)) < 1e-14
    X(end) = [];
    Y(end) = [];
end

% Number of actual panels
Npanel = length(X);

%% ---------------- PANEL GEOMETRY ----------------
x1 = X(:);
y1 = Y(:);
x2 = [X(2:end), X(1)].';
y2 = [Y(2:end), Y(1)].';

dx = x2-x1;
dy = y2-y1;
S  = hypot(dx,dy);
phi = atan2(dy,dx);

% Panel midpoint/control point
xc = 0.5*(x1+x2);
yc = 0.5*(y1+y2);

% For counter-clockwise contour the outward normal is:
% n = [sin(phi), -cos(phi)]
nx = sin(phi);
ny = -cos(phi);

% Tangent direction
tx = cos(phi);
ty = sin(phi);

%% ---------------- FREESTREAM ----------------
alpha = deg2rad(alpha_deg);
Uvec = Uinf*[cos(alpha), sin(alpha)];

%% -------- CONSTANT-SOURCE PANEL INFLUENCE MATRIX --------
% Boundary condition:
%     (V_infinity + V_source) . n = 0
%
% A(i,j) = normal velocity at control point i due to a
%          unit-strength source distributed over panel j.
%
% Diagonal source-panel coefficient = +0.5 for the exterior
% boundary condition with the present counter-clockwise contour.

A = zeros(Npanel,Npanel);
B = zeros(Npanel,1);

for i = 1:Npanel

    % RHS = - freestream normal velocity
    B(i) = -(Uvec(1)*nx(i) + Uvec(2)*ny(i));

    for j = 1:Npanel

        if i == j
            A(i,j) = 0.5;
            continue;
        end

        % Control point i in local coordinates of panel j
        dxp = xc(i)-x1(j);
        dyp = yc(i)-y1(j);

        c = cos(phi(j));
        s = sin(phi(j));

        Xlocal =  dxp*c + dyp*s;
        Ylocal = -dxp*s + dyp*c;

        Sj = S(j);

        r1sq = (Xlocal)^2 + (Ylocal)^2;
        r2sq = (Xlocal-Sj)^2 + (Ylocal)^2;

        % Local velocity induced by a unit constant source panel.
        % Correct signs for the source-panel convention used here.
        uLocal = (1/(4*pi))*log(r1sq/r2sq);
        vLocal = (1/(2*pi))*( ...
            atan2(Ylocal,Xlocal-Sj) - atan2(Ylocal,Xlocal) );

        % Transform to global coordinates
        uGlobal = uLocal*c - vLocal*s;
        vGlobal = uLocal*s + vLocal*c;

        % Dot product with outward normal at control point i
        A(i,j) = uGlobal*nx(i) + vGlobal*ny(i);
    end
end

%% ---------------- SOURCE STRENGTHS ----------------
lambda = A\B;

% Check impermeability residual
normalResidual = A*lambda - B;

fprintf('\n===============================================\n');
fprintf('        TASK 4(c): NACA 0021 SOURCE PANEL\n');
fprintf('===============================================\n');
fprintf('Number of panels              = %d\n',Npanel);
fprintf('NACA thickness ratio          = %.3f\n',t);
fprintf('Angle of attack               = %.2f deg\n',alpha_deg);
fprintf('Maximum |BC residual|         = %.3e\n',max(abs(normalResidual)));
fprintf('Source conservation sum       = %.8e\n',sum(lambda.*S));

%% ---------------- SURFACE VELOCITY & Cp ----------------
Vt = zeros(Npanel,1);

for i = 1:Npanel

    % Freestream tangential component
    Vt(i) = Uvec(1)*tx(i) + Uvec(2)*ty(i);

    for j = 1:Npanel

        if i == j
            % A constant source panel has no finite tangential
            % self-induced velocity at the panel midpoint.
            continue;
        end

        dxp = xc(i)-x1(j);
        dyp = yc(i)-y1(j);

        c = cos(phi(j));
        s = sin(phi(j));

        Xlocal =  dxp*c + dyp*s;
        Ylocal = -dxp*s + dyp*c;

        Sj = S(j);

        r1sq = Xlocal^2 + Ylocal^2;
        r2sq = (Xlocal-Sj)^2 + Ylocal^2;

        uLocal = (1/(4*pi))*log(r1sq/r2sq);
        vLocal = (1/(2*pi))*( ...
            atan2(Ylocal,Xlocal-Sj) - atan2(Ylocal,Xlocal) );

        uGlobal = uLocal*c - vLocal*s;
        vGlobal = uLocal*s + vLocal*c;

        Vt(i) = Vt(i) + lambda(j)*(uGlobal*tx(i) + vGlobal*ty(i));
    end
end

Cp = 1 - (Vt/Uinf).^2;

fprintf('Cp minimum                    = %.8f\n',min(Cp));
fprintf('Cp maximum                    = %.8f\n',max(Cp));

%% ---------------- SPLIT UPPER / LOWER SURFACES ----------------
% With the chosen contour:
% first N/2 panels are upper surface (TE -> LE)
% remaining N/2 panels are lower surface (LE -> TE).

mid = Npanel/2;

xCpUpper = xc(1:mid);
CpUpper  = Cp(1:mid);

xCpLower = xc(mid+1:end);
CpLower  = Cp(mid+1:end);

% Sort by x/c from leading edge to trailing edge for plotting.
[xCpUpper,iu] = sort(xCpUpper,'ascend');
CpUpper = CpUpper(iu);

[xCpLower,il] = sort(xCpLower,'ascend');
CpLower = CpLower(il);

%% ---------------- SYMMETRY CHECK ----------------
% Interpolate both surfaces on a common x-grid.
xCommon = linspace(0.005,0.995,300).';

CpUcommon = interp1(xCpUpper,CpUpper,xCommon,'pchip','extrap');
CpLcommon = interp1(xCpLower,CpLower,xCommon,'pchip','extrap');

symmetryError = max(abs(CpUcommon-CpLcommon));

fprintf('Maximum upper/lower |Cp| diff = %.8e\n',symmetryError);

%% ---------------- XFOIL COMPARISON ----------------
hasXfoil = isfile(xfoilFile);

if hasXfoil

    % Read numeric lines from the XFOIL CpWR output.
    fid = fopen(xfoilFile,'r');
    raw = textscan(fid,'%f %f %f','CommentStyle','#',...
        'Delimiter',{' ','\t'},'MultipleDelimsAsOne',true);
    fclose(fid);

    xfoil_x  = raw{1};
    xfoil_y  = raw{2};
    xfoil_Cp = raw{3};

    valid = isfinite(xfoil_x) & isfinite(xfoil_y) & isfinite(xfoil_Cp);
    xfoil_x  = xfoil_x(valid);
    xfoil_y  = xfoil_y(valid);
    xfoil_Cp = xfoil_Cp(valid);

    % XFOIL usually provides TE -> upper -> LE -> lower -> TE.
    % Determine upper/lower by the sign of y.
    isUpperXF = xfoil_y >= 0;
    isLowerXF = xfoil_y < 0;

    xfu = xfoil_x(isUpperXF);
    Cfu = xfoil_Cp(isUpperXF);

    xfl = xfoil_x(isLowerXF);
    Cfl = xfoil_Cp(isLowerXF);

    [xfu,ku] = sort(xfu,'ascend');
    Cfu = Cfu(ku);

    [xfl,kl] = sort(xfl,'ascend');
    Cfl = Cfl(kl);

    % Interpolate XFOIL onto panel x locations.
    CpU_xf = interp1(xfu,Cfu,xCpUpper,'linear','extrap');
    CpL_xf = interp1(xfl,Cfl,xCpLower,'linear','extrap');

    errU = CpUpper-CpU_xf;
    errL = CpLower-CpL_xf;

    rmseU = sqrt(mean(errU.^2));
    rmseL = sqrt(mean(errL.^2));

    xAll = [xCpUpper(:);xCpLower(:)];
    CpAll = [CpUpper(:);CpLower(:)];
    CpXFall = [CpU_xf(:);CpL_xf(:)];

    rmseAll = sqrt(mean((CpAll-CpXFall).^2));
    maeAll = mean(abs(CpAll-CpXFall));

    fprintf('\n----------- XFOIL COMPARISON -----------\n');
    fprintf('XFOIL file                    = %s\n',xfoilFile);
    fprintf('Upper-surface RMSE            = %.8f\n',rmseU);
    fprintf('Lower-surface RMSE            = %.8f\n',rmseL);
    fprintf('Overall RMSE                  = %.8f\n',rmseAll);
    fprintf('Overall MAE                   = %.8f\n',maeAll);

else
    warning(['XFOIL file not found: ',xfoilFile,...
        '. Panel result will still be plotted.']);
end

%% ---------------- FIGURE: Cp COMPARISON ----------------
fig = figure('Color','w','Position',[100 100 1250 700]);
hold on; box on; grid on;

plot(xCpUpper,CpUpper,'o','MarkerSize',5,...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'DisplayName','64-panel source method');

plot(xCpLower,CpLower,'o','MarkerSize',5,...
    'MarkerFaceColor','k','MarkerEdgeColor','k',...
    'HandleVisibility','off');

if hasXfoil
    plot(xfu,Cfu,'-','LineWidth',1.5,...
        'DisplayName','Inviscid XFOIL');

    plot(xfl,Cfl,'-','LineWidth',1.5,...
        'HandleVisibility','off');
end

set(gca,'YDir','reverse');
xlabel('$x/c$','Interpreter','latex','FontSize',15);
ylabel('$C_p$','Interpreter','latex','FontSize',15);
title('Task 4(c): NACA 0021 Pressure Coefficient at \alpha = 0^\circ',...
    'Interpreter','latex','FontSize',16);

legend('Location','best','FontSize',11);
set(gca,'FontSize',12,'LineWidth',1.0);

xlim([0 1]);

% Save publication-quality figure
exportgraphics(fig,'Task4c_NACA0021_Cp_Panel_vs_XFOIL_alpha0.png',...
    'Resolution',300);

%% ---------------- FIGURE: SOURCE STRENGTH ----------------
fig2 = figure('Color','w','Position',[120 120 1250 650]);
hold on; box on; grid on;

plot(xc(1:mid),lambda(1:mid),'o-','LineWidth',1.2,...
    'DisplayName','Upper surface');
plot(xc(mid+1:end),lambda(mid+1:end),'s-','LineWidth',1.2,...
    'DisplayName','Lower surface');

xlabel('$x/c$','Interpreter','latex','FontSize',15);
ylabel('$\lambda$','Interpreter','latex','FontSize',15);
title('Task 4(c): Constant Source Strength Distribution',...
    'Interpreter','latex','FontSize',16);
legend('Location','best');
set(gca,'FontSize',12,'LineWidth',1.0);
xlim([0 1]);

exportgraphics(fig2,'Task4c_NACA0021_SourceStrength.png',...
    'Resolution',300);

fprintf('\nFigures saved:\n');
fprintf('  Task4c_NACA0021_Cp_Panel_vs_XFOIL_alpha0.png\n');
fprintf('  Task4c_NACA0021_SourceStrength.png\n');
fprintf('===============================================\n');
