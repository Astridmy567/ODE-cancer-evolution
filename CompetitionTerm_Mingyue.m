%% 8-equation, 2-compartment system with motility + one-way mutations
% 4 cell types: H (healthy), X (fast growth), Y (fast motility), Z (fast growth+motility)

%% Time span
tspan = [0 1000];

%% Baseline rates
% growth
rH = 0.1; 
rX = 0.2;      % higher than rH
rY = rH;       % same as H

% motility
fH = 0.01;
fX = fH;       % same as H
fY = 0.05;     % higher than H

% one-way mutation rates
mhx = 0.001;   % H -> X
mhy = 0.001;   % H -> Y
mxz = 0.0001;   % X -> Z
myz = 0.0001;   % Y -> Z

% cross-talk interference (for Z baselines)
ct_r = 0.5;   
ct_f = 0.5;   

rZ = rH + ct_r*(rX - rH);
fZ = fH + ct_f*(fY - fH);

% carrying capacities
KH = 1.0; KX = 1.0; KY = 1.0; KZ = 1.0;

% competition coefficients alpha(i,j): effect of type j on type i
alpha = [1.0  0.5  0.3  0.5;   % effect of H,X,Y,Z on H
         0.3  1.0  0.3  0.5;   % effect of H,X,Y,Z on X
         0.3  0.5  1.0  0.5;   % effect of H,X,Y,Z on Y
         0.3  0.5  0.3  1.0];  % effect of H,X,Y,Z on Z

%% Interclonal cooperativity parameters
% set-up: 1 = linear to population (a*X, b*Y);  
% 2 = saturating (a*X/(k+X), b*Y/(k+Y))
setup = 1;
a = 0.1;            % growth modulation by X
b = 0.1;            % motility modulation by Y
k = 0.1;            

%% Initial conditions
% y = [H0; X0; Y0; Z0; H1; X1; Y1; Z1]
H0_0 = 0.8;  X0_0 = 0.1;  Y0_0 = 0.1;  Z0_0 = 0.0;
H1_0 = 0.0;  X1_0 = 0.0;  Y1_0 = 0.0;  Z1_0 = 0.0;
y0 = [H0_0; X0_0; Y0_0; Z0_0; H1_0; X1_0; Y1_0; Z1_0];

%% ODE solver
[t,Y] = ode45(@(t,y) odesys8_lv(t,y, ...
    rH,rX,rY,rZ, fH,fX,fY,fZ, ...
    mhx,mhy,mxz,myz, ...
    KH,KX,KY,KZ, alpha, ...
    a,b,k, setup), tspan, y0);

H0 = Y(:,1); X0 = Y(:,2); Y0c = Y(:,3); Z0 = Y(:,4);
H1 = Y(:,5); X1 = Y(:,6); Y1c = Y(:,7); Z1 = Y(:,8);

figure; hold on; grid on;
plot(t, H0, '-', 'LineWidth', 2);   plot(t, H1, '--', 'LineWidth', 2);
plot(t, X0, '-', 'LineWidth', 2);   plot(t, X1, '--', 'LineWidth', 2);
plot(t, Y0c,'-', 'LineWidth', 2);   plot(t, Y1c,'--', 'LineWidth', 2);
plot(t, Z0, '-', 'LineWidth', 2);   plot(t, Z1, '--', 'LineWidth', 2);
xlabel('Time'); ylabel('Population');
legend({'H_0','H_1','X_0','X_1','Y_0','Y_1','Z_0','Z_1'}, 'Location','bestoutside');
title('Four types in two compartments (LV competition + motility + mutations + cooperativity)');

%% Final values
C = Y(end,:);  % extract the last row
fprintf('Final populations at t = %.2f\n', t(end));
fprintf('  H0=%.4f  X0=%.4f  Y0=%.4f  Z0=%.4f\n', C(1),C(2),C(3),C(4));
fprintf('  H1=%.4f  X1=%.4f  Y1=%.4f  Z1=%.4f\n', C(5),C(6),C(7),C(8));
fprintf('Totals per type:  H=%.4f  X=%.4f  Y=%.4f  Z=%.4f\n', ...
    C(1)+C(5), C(2)+C(6), C(3)+C(7), C(4)+C(8));
fprintf('Totals per space: N0=%.4f  N1=%.4f\n', sum(C(1:4)), sum(C(5:8)));

%% -------- ODE system with Lotka–Volterra crowding --------
function dydt = odesys8_lv(~, y, ...
    rH,rX,rY,rZ, fH,fX,fY,fZ, ...
    mhx,mhy,mxz,myz, ...
    KH,KX,KY,KZ, alpha, ...
    a,b,k, setup)

    H0 = y(1); X0 = y(2); Y0 = y(3); Z0 = y(4);
    H1 = y(5); X1 = y(6); Y1 = y(7); Z1 = y(8);

    % if statement for setup 1 or 2 (cooperativity)
    if setup == 1
        growthBoostX0 = a*X0;          growthBoostX1 = a*X1;
        motilityBoostY0 = b*Y0;        motilityBoostY1 = b*Y1;
    else
        growthBoostX0 = a*X0/(k+X0);   growthBoostX1 = a*X1/(k+X1);
        motilityBoostY0 = b*Y0/(k+Y0); motilityBoostY1 = b*Y1/(k+Y1);
    end

    % Effective rates per type/space
    rH_p  = rH + growthBoostX0;   rH_pp = rH + growthBoostX1;
    fH_p  = fH + motilityBoostY0; fH_pp = fH + motilityBoostY1;

    rX_p  = rX;                   rX_pp = rX;
    fX_p  = fX + motilityBoostY0; fX_pp = fX + motilityBoostY1;

    rY_p  = rY + growthBoostX0;   rY_pp = rY + growthBoostX1;
    fY_p  = fY;                   fY_pp = fY;

    rZ_p  = rZ;                   rZ_pp = rZ;
    fZ_p  = fZ;                   fZ_pp = fZ;

    % LV competition pressures
    compH0 = (alpha(1,1)*H0 + alpha(1,2)*X0 + alpha(1,3)*Y0 + alpha(1,4)*Z0) / KH;
    compX0 = (alpha(2,1)*H0 + alpha(2,2)*X0 + alpha(2,3)*Y0 + alpha(2,4)*Z0) / KX;
    compY0 = (alpha(3,1)*H0 + alpha(3,2)*X0 + alpha(3,3)*Y0 + alpha(3,4)*Z0) / KY;
    compZ0 = (alpha(4,1)*H0 + alpha(4,2)*X0 + alpha(4,3)*Y0 + alpha(4,4)*Z0) / KZ;

    compH1 = (alpha(1,1)*H1 + alpha(1,2)*X1 + alpha(1,3)*Y1 + alpha(1,4)*Z1) / KH;
    compX1 = (alpha(2,1)*H1 + alpha(2,2)*X1 + alpha(2,3)*Y1 + alpha(2,4)*Z1) / KX;
    compY1 = (alpha(3,1)*H1 + alpha(3,2)*X1 + alpha(3,3)*Y1 + alpha(3,4)*Z1) / KY;
    compZ1 = (alpha(4,1)*H1 + alpha(4,2)*X1 + alpha(4,3)*Y1 + alpha(4,4)*Z1) / KZ;

    % Space 0
    dH0 = rH_p*H0*(1 - compH0) + fH_p*(H1 - H0) - mhx*H0 - mhy*H0;
    dX0 = rX_p*X0*(1 - compX0) + fX_p*(X1 - X0) + mhx*H0 - mxz*X0;
    dY0 = rY_p*Y0*(1 - compY0) + fY_p*(Y1 - Y0) + mhy*H0 - myz*Y0;
    dZ0 = rZ_p*Z0*(1 - compZ0) + fZ_p*(Z1 - Z0) + mxz*X0 + myz*Y0;

    % Space 1
    dH1 = rH_pp*H1*(1 - compH1) + fH_pp*(H0 - H1) - mhx*H1 - mhy*H1;
    dX1 = rX_pp*X1*(1 - compX1) + fX_pp*(X0 - X1) + mhx*H1 - mxz*X1;
    dY1 = rY_pp*Y1*(1 - compY1) + fY_pp*(Y0 - Y1) + mhy*H1 - myz*Y1;
    dZ1 = rZ_pp*Z1*(1 - compZ1) + fZ_pp*(Z0 - Z1) + mxz*X1 + myz*Y1;

    dydt = [dH0; dX0; dY0; dZ0; dH1; dX1; dY1; dZ1];
end