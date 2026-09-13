%% AERO3260 ASSIGNMENT 1 - TASK 4(e)
% FINAL MASTER CODE
%
% NACA 0021, alpha = 10 deg
% 64-panel constant-source + constant-strength vortex panel method
% compared with ACTUAL XFOIL inviscid Cp data.
%
% IMPORTANT:
% The panel-method solver below is deliberately the SAME validated
% formulation used for Task 4(d). It is not a new/revised solver.
%
% XFOIL data file:
%       xfoil_cp10.txt
% The script first checks the script folder and Current Folder.
% If the file is elsewhere, a file-selection window is opened automatically.
%
% XFOIL settings used to generate that file:
%       NACA 0021
%       PANE
%       OPER
%       (inviscid .OPERi mode)
%       ALFA 10
%       CPWR xfoil_cp10.txt
%
% XFOIL run reported:
%       alpha = 10 deg
%       CL    = 1.2870
%       CM    = -0.0297
%       Re    = 0 (inviscid)
%
% This script:
%   1. Rebuilds the validated 64-panel vortex-panel solution.
%   2. Checks source conservation, Kutta residual and conditioning.
%   3. Reads the actual XFOIL Cp file.
%   4. Plots panel and XFOIL Cp on one clean figure.
%   5. Calculates Cp comparison errors by surface.
%   6. Plots the geometry comparison.
%   7. Saves report-ready figures and CSV data.
%
% -----------------------------------------------------------------

clear;
clc;
close all;

%% ================================================================
% 1. INPUTS
% ================================================================

N       = 64;
t       = 0.21;
c       = 1.0;
Uinf    = 1.0;
alphaDeg = 10.0;
alpha   = deg2rad(alphaDeg);

xfoilFile = 'xfoil_cp10.txt';

% Locate the XFOIL data robustly. This avoids requiring the MATLAB
% Current Folder to be the same as the folder containing this script.
scriptFolder = fileparts(mfilename('fullpath'));
candidates = {fullfile(scriptFolder,xfoilFile), fullfile(pwd,xfoilFile)};
foundFile = '';
for k = 1:numel(candidates)
    if isfile(candidates{k})
        foundFile = candidates{k};
        break;
    end
end

if isempty(foundFile)
    [f,p] = uigetfile({'*.txt','XFOIL text files (*.txt)'}, ...
        'Select the XFOIL Cp file (e.g. xfoil_cp10.txt)');
    if isequal(f,0)
        error('No XFOIL Cp file was selected.');
    end
    foundFile = fullfile(p,f);
end

xfoilFile = foundFile;

% XFOIL integrated result from the actual converged run shown in XFOIL
CL_XFOIL = 1.2870;
CM_XFOIL = -0.0297;

fprintf('\n============================================================\n');
fprintf('AERO3260 TASK 4(e) - NACA 0021 PANEL METHOD vs XFOIL\n');
fprintf('============================================================\n');
fprintf('NACA 0021\n');
fprintf('Panels                  = %d\n',N);
fprintf('Angle of attack         = %.2f deg\n',alphaDeg);
fprintf('Freestream velocity     = %.4f\n',Uinf);
fprintf('XFOIL file              = %s\n',xfoilFile);

%% ================================================================
% 2. GENERATE NACA 0021 GEOMETRY
% ================================================================

% Cosine-spaced x coordinates, LE -> TE
xHalf = (1-cos(linspace(0,pi,N/2+1)))/2;

% NACA 0021 thickness distribution
yt = 5*t*( ...
      0.2969*sqrt(xHalf) ...
    - 0.1260*xHalf ...
    - 0.3516*xHalf.^2 ...
    + 0.2843*xHalf.^3 ...
    - 0.1015*xHalf.^4 );

% Close TE for consistency with the validated Task 4(c)/(d) MATLAB
% geometry.
yt(end) = 0;

% Boundary ordering:
% upper TE -> LE -> lower TE
x = [fliplr(xHalf), xHalf(2:end)];
y = [fliplr(yt),    -yt(2:end)];

x = x(:);
y = y(:);

if numel(x) ~= N+1
    error('Geometry error: expected N+1 nodes.');
end

%% ================================================================
% 3. PANEL GEOMETRY
% ================================================================

x1 = x(1:N);
y1 = y(1:N);
x2 = x(2:N+1);
y2 = y(2:N+1);

dx = x2-x1;
dy = y2-y1;

S   = hypot(dx,dy);
phi = atan2(dy,dx);

xc = 0.5*(x1+x2);
yc = 0.5*(y1+y2);

% CLOCKWISE panel ordering.
% These are the inward normals used in the validated Task 4(d)
% formulation.
nx = -sin(phi);
ny =  cos(phi);

% Panel tangents
tx = cos(phi);
ty = sin(phi);

%% ================================================================
% 4. SOURCE + VORTEX INFLUENCE MATRICES
% ================================================================

As = zeros(N,N);      % source -> normal
Ts = zeros(N,N);      % source -> tangential

Av = zeros(N,1);      % constant vortex -> normal
Tv = zeros(N,1);      % constant vortex -> tangential

for i = 1:N
    for j = 1:N

        % Control point i in panel-j local coordinates
        X = (xc(i)-x1(j))*cos(phi(j)) ...
          + (yc(i)-y1(j))*sin(phi(j));

        Y = -(xc(i)-x1(j))*sin(phi(j)) ...
          + (yc(i)-y1(j))*cos(phi(j));

        if i == j

            % Validated self influence for the chosen clockwise
            % boundary/normal convention.
            uLocal = 0.0;
            vLocal = -0.5;

        else

            r1sq = X^2 + Y^2;
            r2sq = (X-S(j))^2 + Y^2;

            uLocal = (1/(4*pi))*log(r1sq/r2sq);

            dTheta = atan2(Y,X-S(j)) - atan2(Y,X);

            % Prevent atan2 branch jumps
            dTheta = atan2(sin(dTheta),cos(dTheta));

            vLocal = (1/(2*pi))*dTheta;
        end

        % Source-panel velocity in global coordinates
        uS = uLocal*cos(phi(j)) - vLocal*sin(phi(j));
        vS = uLocal*sin(phi(j)) + vLocal*cos(phi(j));

        % Source normal/tangential influence
        As(i,j) = uS*nx(i) + vS*ny(i);
        Ts(i,j) = uS*tx(i) + vS*ty(i);

        % Corresponding vortex-panel velocity
        uV = -vS;
        vV =  uS;

        Av(i) = Av(i) + uV*nx(i) + vV*ny(i);
        Tv(i) = Tv(i) + uV*tx(i) + vV*ty(i);

    end
end

%% ================================================================
% 5. FREESTREAM
% ================================================================

uInf = Uinf*cos(alpha);
vInf = Uinf*sin(alpha);

VnInf = uInf*nx + vInf*ny;
VtInf = uInf*tx + vInf*ty;

%% ================================================================
% 6. 65 x 65 SYSTEM
% ================================================================
% Unknowns:
% lambda(1)...lambda(64), gamma
%
% 64 no-penetration equations
% 1 Kutta equation

A = zeros(N+1,N+1);
b = zeros(N+1,1);

A(1:N,1:N) = As;
A(1:N,N+1) = Av;
b(1:N)     = -VnInf;

% Kutta condition at the trailing edge:
% Vt(panel 1) + Vt(panel N) = 0
A(N+1,1:N) = Ts(1,:) + Ts(N,:);
A(N+1,N+1) = Tv(1) + Tv(N);

b(N+1) = -(VtInf(1) + VtInf(N));

%% ================================================================
% 7. SOLVE
% ================================================================

q = A\b;

lambda = q(1:N);
gamma  = q(N+1);

%% ================================================================
% 8. SURFACE TANGENTIAL VELOCITY AND Cp
% ================================================================

Vt = VtInf + Ts*lambda + Tv*gamma;

CpPanel = 1 - (Vt/Uinf).^2;

%% ================================================================
% 9. NUMERICAL VALIDATION
% ================================================================

% Source conservation
sourceResidual = sum(lambda.*S);

% Kutta residual
kuttaResidual = Vt(1) + Vt(N);

% Matrix conditioning
conditionNumber = cond(A);

CpMinPanel = min(CpPanel);
CpMaxPanel = max(CpPanel);

%% ================================================================
% 10. AERODYNAMIC FORCE COEFFICIENTS
% ================================================================
% The stored nx,ny are inward normals.
% Pressure force = Cp * [nx,ny] * ds.
%
% First calculate body-axis coefficients:
Cx = sum(CpPanel .* nx .* S/c);
Cy = sum(CpPanel .* ny .* S/c);

% Convert to lift/drag relative to freestream.
CL_panel = -Cx*sin(alpha) + Cy*cos(alpha);
CD_panel =  Cx*cos(alpha) + Cy*sin(alpha);

%% ================================================================
% 11. SPLIT PANEL Cp INTO UPPER AND LOWER SURFACES
% ================================================================

% With this boundary ordering:
% panels 1:N/2   = upper surface, TE -> LE
% panels N/2+1:N = lower surface, LE -> TE

idxU = 1:N/2;
idxL = N/2+1:N;

xUpperPanel  = xc(idxU);
CpUpperPanel = CpPanel(idxU);

xLowerPanel  = xc(idxL);
CpLowerPanel = CpPanel(idxL);

% Sort both surfaces from LE -> TE
[xUpperPanel,ordU] = sort(xUpperPanel);
CpUpperPanel = CpUpperPanel(ordU);

[xLowerPanel,ordL] = sort(xLowerPanel);
CpLowerPanel = CpLowerPanel(ordL);

%% ================================================================
% 12. READ ACTUAL XFOIL FILE
% ================================================================

if ~isfile(xfoilFile)
    error('XFOIL Cp file could not be located.');
end

% The actual XFOIL file has:
% line 1: airfoil name
% line 2: alpha / Re / etc.
% line 3: column headings
% data starts after these lines.
xfoilData = readmatrix(xfoilFile,...
    'FileType','text',...
    'NumHeaderLines',3);

if size(xfoilData,2) < 3
    error('XFOIL file must contain at least x, y and Cp columns.');
end

xXfoil  = xfoilData(:,1);
yXfoil  = xfoilData(:,2);
CpXfoil = xfoilData(:,3);

% Remove non-numeric rows if present
valid = isfinite(xXfoil) & isfinite(yXfoil) & isfinite(CpXfoil);

xXfoil  = xXfoil(valid);
yXfoil  = yXfoil(valid);
CpXfoil = CpXfoil(valid);

if isempty(xXfoil)
    error('No valid XFOIL x-y-Cp data were found.');
end

fprintf('\nXFOIL points read = %d\n',length(xXfoil));

%% ================================================================
% 13. SPLIT XFOIL UPPER / LOWER SURFACES
% ================================================================

% Use y sign, while assigning exact LE/TE points safely.
upperXF = yXfoil >= 0;
lowerXF = yXfoil <  0;

xUpperXfoil  = xXfoil(upperXF);
CpUpperXfoil = CpXfoil(upperXF);

xLowerXfoil  = xXfoil(lowerXF);
CpLowerXfoil = CpXfoil(lowerXF);

% Sort LE -> TE
[xUpperXfoil,ordU] = sort(xUpperXfoil);
CpUpperXfoil = CpUpperXfoil(ordU);

[xLowerXfoil,ordL] = sort(xLowerXfoil);
CpLowerXfoil = CpLowerXfoil(ordL);

%% ================================================================
% 14. XFOIL DATA CHECK
% ================================================================

CpMinXfoil = min(CpXfoil);
CpMaxXfoil = max(CpXfoil);

%% ================================================================
% 15. Cp ERROR ANALYSIS
% ================================================================
% Compare XFOIL and panel Cp at the panel control-point x locations.
% Interpolation is performed separately on upper and lower surfaces.

% Remove duplicate x values before interpolation if necessary.
[xUpperXfoil,uniqueU] = unique(xUpperXfoil,'stable');
CpUpperXfoil = CpUpperXfoil(uniqueU);

[xLowerXfoil,uniqueL] = unique(xLowerXfoil,'stable');
CpLowerXfoil = CpLowerXfoil(uniqueL);

% Interpolate XFOIL Cp at panel control points
CpUpperXF_interp = interp1( ...
    xUpperXfoil,CpUpperXfoil,...
    xUpperPanel,'linear','extrap');

CpLowerXF_interp = interp1( ...
    xLowerXfoil,CpLowerXfoil,...
    xLowerPanel,'linear','extrap');

% Errors
errUpper = CpUpperPanel - CpUpperXF_interp;
errLower = CpLowerPanel - CpLowerXF_interp;

RMSE_upper = sqrt(mean(errUpper.^2));
RMSE_lower = sqrt(mean(errLower.^2));

MAE_upper = mean(abs(errUpper));
MAE_lower = mean(abs(errLower));

RMSE_all = sqrt( ...
    mean([errUpper(:);errLower(:)].^2));

MAE_all = mean(abs([errUpper(:);errLower(:)]));

%% ================================================================
% 16. CL COMPARISON
% ================================================================

CL_difference = CL_panel - CL_XFOIL;
CL_percentDifference = ...
    abs(CL_difference)/abs(CL_XFOIL)*100;

%% ================================================================
% 17. PRINT ALL NUMERICAL CHECKS
% ================================================================

fprintf('\n============================================================\n');
fprintf('VALIDATION OF 64-PANEL VORTEX-PANEL SOLUTION\n');
fprintf('============================================================\n');

fprintf('Gamma (constant vortex strength) = %.10f\n',gamma);
fprintf('Source conservation residual     = %.8e\n',sourceResidual);
fprintf('Kutta residual                   = %.8e\n',kuttaResidual);
fprintf('Matrix condition number          = %.6e\n',conditionNumber);

fprintf('\nPanel-method Cp:\n');
fprintf('Cp minimum                       = %.8f\n',CpMinPanel);
fprintf('Cp maximum                       = %.8f\n',CpMaxPanel);

fprintf('\nPanel-method forces:\n');
fprintf('CL                               = %.8f\n',CL_panel);
fprintf('CD                               = %.8f\n',CD_panel);

fprintf('\nXFOIL data:\n');
fprintf('Number of XFOIL points           = %d\n',length(xXfoil));
fprintf('Cp minimum                       = %.8f\n',CpMinXfoil);
fprintf('Cp maximum                       = %.8f\n',CpMaxXfoil);
fprintf('CL from XFOIL run                = %.8f\n',CL_XFOIL);
fprintf('CM from XFOIL run                = %.8f\n',CM_XFOIL);

fprintf('\nCp comparison:\n');
fprintf('Upper-surface RMSE               = %.8f\n',RMSE_upper);
fprintf('Lower-surface RMSE               = %.8f\n',RMSE_lower);
fprintf('Overall Cp RMSE                  = %.8f\n',RMSE_all);
fprintf('Upper-surface MAE                = %.8f\n',MAE_upper);
fprintf('Lower-surface MAE                = %.8f\n',MAE_lower);
fprintf('Overall Cp MAE                   = %.8f\n',MAE_all);

fprintf('\nCL comparison:\n');
fprintf('CL panel                         = %.8f\n',CL_panel);
fprintf('CL XFOIL                         = %.8f\n',CL_XFOIL);
fprintf('CL difference                    = %.8f\n',CL_difference);
fprintf('CL percentage difference         = %.4f %%\n',CL_percentDifference);

fprintf('============================================================\n');

%% ================================================================
% 18. FIGURE 1 - FINAL Cp COMPARISON
% ================================================================

fig1 = figure( ...
    'Color','w',...
    'Position',[100 100 1200 760],...
    'Name','Task 4(e) - Cp Comparison');

hold on;
box on;
grid on;

% XFOIL first: smooth continuous curves
plot(xUpperXfoil,CpUpperXfoil,...
    'b-',...
    'LineWidth',2.0,...
    'DisplayName','XFOIL - upper surface');

plot(xLowerXfoil,CpLowerXfoil,...
    'r-',...
    'LineWidth',2.0,...
    'DisplayName','XFOIL - lower surface');

% Panel method: control-point markers
plot(xUpperPanel,CpUpperPanel,...
    'ko-',...
    'LineWidth',1.1,...
    'MarkerSize',4,...
    'MarkerFaceColor','k',...
    'DisplayName','64-panel method - upper surface');

plot(xLowerPanel,CpLowerPanel,...
    'ks-',...
    'LineWidth',1.1,...
    'MarkerSize',4,...
    'MarkerFaceColor','w',...
    'DisplayName','64-panel method - lower surface');

% Conventional Cp presentation
set(gca,'YDir','reverse');

xlabel('x/c','FontSize',14);
ylabel('C_p','FontSize',14);

title('NACA 0021: 64-Panel Method vs XFOIL, \alpha = 10^\circ',...
    'FontSize',16);

legend('Location','best','FontSize',11);

xlim([0 1]);

% Let MATLAB choose a sensible y range from both datasets
allCp = [CpPanel(:);CpXfoil(:)];
cpLo = min(allCp);
cpHi = max(allCp);
margin = 0.08*(cpHi-cpLo);
ylim([cpLo-margin cpHi+margin]);

set(gca,'FontSize',12,'LineWidth',1);

%% ================================================================
% 19. SAVE FINAL Cp FIGURE
% ================================================================

exportgraphics(fig1,...
    'Task4e_NACA0021_Cp_Panel_vs_XFOIL.png',...
    'Resolution',300);

%% ================================================================
% 20. FIGURE 2 - GEOMETRY CHECK
% ================================================================

fig2 = figure( ...
    'Color','w',...
    'Position',[120 120 1100 650],...
    'Name','Task 4(e) - Geometry Check');

hold on;
box on;
grid on;

plot(x,y,...
    'k-',...
    'LineWidth',1.7,...
    'DisplayName','NACA 0021 panel geometry');

plot(xc,yc,...
    'ko',...
    'MarkerSize',4,...
    'MarkerFaceColor','w',...
    'DisplayName','64 panel control points');

plot(xXfoil,yXfoil,...
    'b.',...
    'MarkerSize',6,...
    'DisplayName','XFOIL points');

axis equal;
xlim([-0.03 1.03]);
ylim([-0.14 0.14]);

xlabel('x/c','FontSize',14);
ylabel('y/c','FontSize',14);

title('NACA 0021 Geometry: Panel Control Points and XFOIL Points',...
    'FontSize',16);

legend('Location','best','FontSize',11);

set(gca,'FontSize',12,'LineWidth',1);

%% ================================================================
% 21. SAVE GEOMETRY FIGURE
% ================================================================

exportgraphics(fig2,...
    'Task4e_NACA0021_Geometry_Check.png',...
    'Resolution',300);

%% ================================================================
% 22. EXPORT PANEL DATA
% ================================================================

panelTable = table( ...
    xc/c,...
    yc/c,...
    CpPanel,...
    Vt/Uinf,...
    lambda,...
    S/c,...
    'VariableNames',{...
    'x_c',...
    'y_c',...
    'Cp_panel',...
    'Vt_over_Uinf',...
    'source_strength',...
    'panel_length_c'});

writetable(panelTable,...
    'Task4e_Panel_Method_Data.csv');

%% ================================================================
% 23. EXPORT XFOIL DATA
% ================================================================

xfoilTable = table( ...
    xXfoil,...
    yXfoil,...
    CpXfoil,...
    'VariableNames',{...
    'x_c',...
    'y_c',...
    'Cp_XFOIL'});

writetable(xfoilTable,...
    'Task4e_XFOIL_Data.csv');

%% ================================================================
% 24. EXPORT COMPARISON ERROR DATA
% ================================================================

comparisonTable = table( ...
    xUpperPanel,...
    CpUpperPanel,...
    CpUpperXF_interp,...
    errUpper,...
    'VariableNames',{...
    'x_c_upper',...
    'Cp_panel_upper',...
    'Cp_XFOIL_upper',...
    'Cp_error_upper'});

writetable(comparisonTable,...
    'Task4e_UpperSurface_Comparison.csv');

comparisonTable2 = table( ...
    xLowerPanel,...
    CpLowerPanel,...
    CpLowerXF_interp,...
    errLower,...
    'VariableNames',{...
    'x_c_lower',...
    'Cp_panel_lower',...
    'Cp_XFOIL_lower',...
    'Cp_error_lower'});

writetable(comparisonTable2,...
    'Task4e_LowerSurface_Comparison.csv');

%% ================================================================
% 25. FINAL STATUS
% ================================================================

fprintf('\n============================================================\n');
fprintf('TASK 4(e) PROCESS FINISHED\n');
fprintf('============================================================\n');

fprintf('Saved figure:\n');
fprintf('  Task4e_NACA0021_Cp_Panel_vs_XFOIL.png\n');

fprintf('Saved figure:\n');
fprintf('  Task4e_NACA0021_Geometry_Check.png\n');

fprintf('Saved data files:\n');
fprintf('  Task4e_Panel_Method_Data.csv\n');
fprintf('  Task4e_XFOIL_Data.csv\n');
fprintf('  Task4e_UpperSurface_Comparison.csv\n');
fprintf('  Task4e_LowerSurface_Comparison.csv\n');

fprintf('\nCheck the Kutta residual above.\n');
fprintf('A value close to zero confirms the Kutta condition.\n');
fprintf('============================================================\n');
