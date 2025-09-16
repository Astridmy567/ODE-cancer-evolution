%% 8-equation, 2-compartment system with motility + one-way mutations
% No LV alpha-competition; uses -mu*T*state crowding

clear; clc;

%% Time span
tspan = [0 1000];

%% Baseline rates
% growth
rH = 0.1;
rX = 0.2;      % higher than rH
rY = rH;       % same as H

% motility
fH = 0;        % 0.01 typical; set 0 if you want to isolate effects
fX = fH;       % same as H
fY = 0;        % 0.05 typical; set 0 here per your snippet

% one-way mutation rates
mhx = 0.001;   % H -> X
mhy = 0.001;   % H -> Y
mxz = 0.001;   % X -> Z
myz = 0.001;   % Y -> Z

% cross-talk interpolation for Z baselines
ct_r = 0.5;      % 0 -> rZ=rH; 1 -> rZ=rX
ct_f = 0.5      % 0 -> fZ=fH; 1 -> fZ=fY

rZ = rH + ct_r*(rX - rH);
fZ = fH + ct_f*(fY - fH);

% crowding strength for -mu*T*state
mu = 0.2;  

%% Interclonal cooperativity parameters (currently OFF with a=b=0)
% setup: 1 = linear (a*X, b*Y); 2 = saturating (a*X/(k+X), b*Y/(k+Y))
setup = 2;
a = 0.2;         % growth modulation by X
b = 0.3;         % motility modulation by Y
k = 0.1;

%% Initial conditions
% y = [H0; X0; Y0; Z0; H1; X1; Y1; Z1]
H0_0 = 0.7;  X0_0 = 0.1;  Y0_0 = 0.1;  Z0_0 = 0.0;
H1_0 = 0.0;  X1_0 = 0.0;  Y1_0 = 0.0;  Z1_0 = 0.0;
y0 = [H0_0; X0_0; Y0_0; Z0_0; H1_0; X1_0; Y1_0; Z1_0];

%% ODE solver
[t,Y] = ode45(@(t,y) odesys8_bd(t,y, ...
    rH,rX,rY,rZ, fH,fX,fY,fZ, ...
    mhx,mhy,mxz,myz, ...
    a,b,k, setup, mu), ...
    tspan, y0);

%% Extract
H0 = Y(:,1); X0 = Y(:,2); Y0 = Y(:,3); Z0 = Y(:,4);
H1 = Y(:,5); X1 = Y(:,6); Y1 = Y(:,7); Z1 = Y(:,8);

%% Plot
figure; hold on; grid on;
plot(t, H0, '-',  'LineWidth', 2);  plot(t, H1, '--', 'LineWidth', 2);
plot(t, X0, '-',  'LineWidth', 2);  plot(t, X1, '--', 'LineWidth', 2);
plot(t, Y0, '-',  'LineWidth', 2);  plot(t, Y1, '--', 'LineWidth', 2);
plot(t, Z0, '-',  'LineWidth', 2);  plot(t, Z1, '--', 'LineWidth', 2);
xlabel('Time'); ylabel('Population');
legend({'H_0','H_1','X_0','X_1','Y_0','Y_1','Z_0','Z_1'}, 'Location','bestoutside');
title('Four types in two compartments (-\mu T crowding + motility + mutations)');

%% Final values
C = Y(end,:);
fprintf('Final populations at t = %.2f\n', t(end));
fprintf('  H0=%.4f  X0=%.4f  Y0=%.4f  Z0=%.4f\n', C(1),C(2),C(3),C(4));
fprintf('  H1=%.4f  X1=%.4f  Y1=%.4f  Z1=%.4f\n', C(5),C(6),C(7),C(8));
fprintf('Totals per type:  H=%.4f  X=%.4f  Y=%.4f  Z=%.4f\n', ...
    C(1)+C(5), C(2)+C(6), C(3)+C(7), C(4)+C(8));
fprintf('Totals per space: N0=%.4f  N1=%.4f\n', sum(C(1:4)), sum(C(5:8)));

%% ODE system: birth-death process
% y = [H0; X0; Y0; Z0; H1; X1; Y1; Z1]
function dydt = odesys8_bd(~, y, ...
    rH,rX,rY,rZ, fH,fX,fY,fZ, ...
    mhx,mhy,mxz,myz, ...
    a,b,k, setup, mu)

    H0 = y(1); X0 = y(2); Y0 = y(3); Z0 = y(4);
    H1 = y(5); X1 = y(6); Y1 = y(7); Z1 = y(8);

    % Interclonal cooperativity
    if setup == 1
        % linear
        growthBoostX0   = a*X0;           growthBoostX1   = a*X1;
        motilityBoostY0 = b*Y0;           motilityBoostY1 = b*Y1;
    else
        % saturating
        growthBoostX0   = a*X0/(k+X0);    growthBoostX1   = a*X1/(k+X1);
        motilityBoostY0 = b*Y0/(k+Y0);    motilityBoostY1 = b*Y1/(k+Y1);
    end

    % Effective rates
    rH_p  = rH + growthBoostX0;   rH_pp = rH + growthBoostX1;
    fH_p  = fH + motilityBoostY0; fH_pp = fH + motilityBoostY1;

    rX_p  = rX;                   rX_pp = rX;
    fX_p  = fX + motilityBoostY0; fX_pp = fX + motilityBoostY1;

    rY_p  = rY + growthBoostX0;   rY_pp = rY + growthBoostX1;
    fY_p  = fY;                   fY_pp = fY;

    % Z not cooperatively boosted
    rZ_p  = rZ;                   rZ_pp = rZ;
    fZ_p  = fZ;                   fZ_pp = fZ;
    
    

    % ODEs
    % Space 0
    T0 = H0 + X0 + Y0 + Z0;
    dH0 = rH_p*H0 - mu*T0*H0 + fH_p*(H1 - H0) - mhx*(H0 - 0*X0) - mhy*(H0 - 0*Y0);
    dX0 = rX_p*X0 - mu*T0*X0 + fX_p*(X1 - X0) + mhx*(H0 - 0*X0) - mxz*(X0 - 0*Z0);
    dY0 = rY_p*Y0 - mu*T0*Y0 + fY_p*(Y1 - Y0) + mhy*(H0 - 0*Y0) - myz*(Y0 - 0*Z0);
    dZ0 = rZ_p*Z0 - mu*T0*Z0 + fZ_p*(Z1 - Z0) + mxz*(X0 - 0*Z0) + myz*(Y0 - 0*Z0);

    % Space 1
    T1 = H1 + X1 + Y1 + Z1;
    dH1 = rH_pp*H1 - mu*T1*H1 + fH_pp*(H0 - H1) - mhx*(H1 - 0*X1) - mhy*(H1 - 0*Y1);
    dX1 = rX_pp*X1 - mu*T1*X1 + fX_pp*(X0 - X1) + mhx*(H1 - 0*X1) - mxz*(X1 - 0*Z1);
    dY1 = rY_pp*Y1 - mu*T1*Y1 + fY_pp*(Y0 - Y1) + mhy*(H1 - 0*Y1) - myz*(Y1 - 0*Z1);
    dZ1 = rZ_pp*Z1 - mu*T1*Z1 + fZ_pp*(Z0 - Z1) + mxz*(X1 - 0*Z1) + myz*(Y1 - 0*Z1);

    dydt = [dH0; dX0; dY0; dZ0; dH1; dX1; dY1; dZ1];
end