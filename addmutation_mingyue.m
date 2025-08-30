%% 8-equation, 2-compartment system with motility + one-way mutations
% 4 cell types: H (healthy), X (fast growth), Y (fast motility), Z (fast growth+motility)

%% Time span
tspan = [0 200];

%% Rates
% growth
rH = 0.1; 
rX = 0.2;      % higher than rH
rY = rH;       % same as H
rZ = rX;       % same as X

% motility
fH = 0.01;
fX = fH;       % same as H
fY = 0.05;     % higher than H
fZ = fY;       % same as Y

% one-way mutation rates
mhx = 0.001;   % H -> X
mhy = 0.001;   % H -> Y
mxz = 0.001;   % X -> Z
myz = 0.001;   % Y -> Z

% capacity
K = 1.0;

%% Initial conditions
H0_0 = 0.8;  X0_0 = 0.1;  Y0_0 = 0.1;  Z0_0 = 0.0;
H1_0 = 0.0;  X1_0 = 0.0;  Y1_0 = 0.0;  Z1_0 = 0.0;
y0 = [H0_0; X0_0; Y0_0; Z0_0; H1_0; X1_0; Y1_0; Z1_0];

%% Ode solver
[t,Y] = ode45(@(t,y) odesys8(t,y, ...
    rH,rX,rY,rZ, fH,fX,fY,fZ, mhx,mhy,mxz,myz, K), tspan, y0);

%% Extract solutions
H0 = Y(:,1); X0 = Y(:,2); Y0c = Y(:,3); Z0 = Y(:,4);
H1 = Y(:,5); X1 = Y(:,6); Y1c = Y(:,7); Z1 = Y(:,8);

%% Plot
figure; hold on; grid on;
plot(t, H0, '-', 'LineWidth', 2);   plot(t, H1, '--', 'LineWidth', 2);
plot(t, X0, '-', 'LineWidth', 2);   plot(t, X1, '--', 'LineWidth', 2);
plot(t, Y0c,'-', 'LineWidth', 2);   plot(t, Y1c,'--', 'LineWidth', 2);
plot(t, Z0, '-', 'LineWidth', 2);   plot(t, Z1, '--', 'LineWidth', 2);
xlabel('Time'); ylabel('Population');
legend({'H_0','H_1','X_0','X_1','Y_0','Y_1','Z_0','Z_1'}, 'Location','bestoutside');
title('Four cell types in two compartments (logistic growth + motility + mutations)');

%% final values
C = Y(end,:);  % last row
fprintf('Final populations at t = %.2f\n', t(end));
fprintf('  H0=%.4f  X0=%.4f  Y0=%.4f  Z0=%.4f\n', C(1),C(2),C(3),C(4));
fprintf('  H1=%.4f  X1=%.4f  Y1=%.4f  Z1=%.4f\n', C(5),C(6),C(7),C(8));
fprintf('Totals per type:  H=%.4f  X=%.4f  Y=%.4f  Z=%.4f\n', ...
    C(1)+C(5), C(2)+C(6), C(3)+C(7), C(4)+C(8));
fprintf('Totals per compartment: N0=%.4f  N1=%.4f\n', sum(C(1:4)), sum(C(5:8)));

%% -------- ODE system --------
% y = [H0; X0; Y0; Z0; H1; X1; Y1; Z1]
function dydt = odesys8(~, y, rH, rX, rY, rZ, fH, fX, fY, fZ, mhx, mhy, mxz, myz, K)

    H0 = y(1); X0 = y(2); Y0 = y(3); Z0 = y(4);
    H1 = y(5); X1 = y(6); Y1 = y(7); Z1 = y(8);

    % Total population in each space
    N0 = H0 + X0 + Y0 + Z0;
    N1 = H1 + X1 + Y1 + Z1;

    % Space 0
    dH0 = rH*H0*(1 - N0/K) + fH*(H1 - H0) - mhx*H0 - mhy*H0;
    dX0 = rX*X0*(1 - N0/K) + fX*(X1 - X0) + mhx*H0 - mxz*X0;
    dY0 = rY*Y0*(1 - N0/K) + fY*(Y1 - Y0) + mhy*H0 - myz*Y0;
    dZ0 = rZ*Z0*(1 - N0/K) + fZ*(Z1 - Z0) + mxz*X0 + myz*Y0;

    % Space 1
    dH1 = rH*H1*(1 - N1/K) + fH*(H0 - H1) - mhx*H1 - mhy*H1;
    dX1 = rX*X1*(1 - N1/K) + fX*(X0 - X1) + mhx*H1 - mxz*X1;
    dY1 = rY*Y1*(1 - N1/K) + fY*(Y0 - Y1) + mhy*H1 - myz*Y1;
    dZ1 = rZ*Z1*(1 - N1/K) + fZ*(Z0 - Z1) + mxz*X1 + myz*Y1;

    dydt = [dH0; dX0; dY0; dZ0; dH1; dX1; dY1; dZ1];
end