% base model
%% single logistic equation - one population
tspan = [0 100];
y0 = 0.1;
[t,y] = ode45(@(t,y) 0.1*y*(1-y), tspan, y0);

figure
plot(t, y, '-o')
xlabel('Time t')
ylabel('Population H(t)')
title('Logistic Growth: r = 0.1, K = 1')
grid on


%% logistic system - two competing populations

% H(0) = 0.1, C(0) = 0.05
y0 = [0.1; 0.05];

f = @(t, y) [ 0.1*y(1)*(1 - (y(1) + y(2))) ;   % dH/dt
              0.5*y(2)*(1 - (y(1) + y(2)))];   % dC/dt

[t, Y] = ode45(f, tspan, y0);


H = Y(:,1); % extract solutions for H from the first column
C = Y(:,2);

figure
plot(t, H, 'b-', 'LineWidth', 2); hold on
plot(t, C, 'r-', 'LineWidth', 2);
xlabel('Time t')
ylabel('Population')
legend('H(t) - Healthy Cells', 'C(t) - Cancer Cells')
title('Logistic Growth System: r_H = 0.1, r_C = 0.5, K = 1')
grid on

%% add motility - four-equation system

% healthy cell growth rates: r_H = 0.1
% cancer cell growth rate: r_C = 0.5
% healthy cell flow rate: f_H = 0.01
% cancer cell flow rate: f_C = 0.05

% H0: y(1), C0: y(2), H1: y(3), C1: y(4); K=0
y0 = [0.1; 0.05; 0.2; 0.06];

f = @(t, y) [ 0.1*y(1)*(1 - (y(1) + y(2))) + 0.01*(y(3)-y(1)) ;   % dH0/dt
              0.5*y(2)*(1 - (y(1) + y(2))) + 0.05*(y(4)-y(2)) ;   % dC0/dt
              0.1*y(3)*(1 - (y(3) + y(4))) + 0.01*(y(1)-y(3)) ;   % dH1/dt
              0.5*y(4)*(1 - (y(3) + y(4))) + 0.05*(y(2)-y(4)) ;]; % dC1/dt 

[t, Y] = ode45(f, tspan, y0);


H0 = Y(:,1); % extract solutions for H from the first column
C0 = Y(:,2);
H1 = Y(:,3);
C1 = Y(:,4);

figure
plot(t, H0, 'b-', 'LineWidth', 2); hold on
plot(t, C0, 'r-', 'LineWidth', 2); hold on
plot(t, H1, 'b-', 'LineWidth', 4); hold on
plot(t, C1, 'r-', 'LineWidth', 4); 
xlabel('Time t')
ylabel('Population')
legend('H_i(t) - Healthy Cells', 'C_i(t) - Cancer Cells')
title('Logistic Growth System')
grid on
