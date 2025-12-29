clc;
clear all;
close all;


% Data
x = [-0.15 -0.86 -1.09]';
y1 = [-0.31 -0.4 -0.12]';%for HCOOH
y2 = [-0.31 -0.38 -0.42]';%for HCHO
y3= [-0.31 -0.38 -0.40]';%for CH3OH
y4= [-0.31 -0.38 -0.28]';%for CH4


% Least-squares linear fitting
p1 = polyfit(x, y1, 1);   % p(1)=slope, p(2)=intercept
p2 = polyfit(x, y2, 1);
p3 = polyfit(x, y3, 1);
p4 = polyfit(x, y4, 1);

% Extract parameters
m1 = -0.1243;
c1 = -0.6;

m2 = p2(1);
c2 = p2(2);

m3 = p3(1);
c3 = p3(2);

m4 = p4(1);
c4 = p4(2);

% Fitted values
y_fit1 = polyval(p1, x);
y_fit2 = polyval(p2, x);
y_fit3 = polyval(p3, x);
y_fit4 = polyval(p4, x);

% % Calculate mean squared error (MSE)
% MSE = mean((y - y_fit).^2);


% ----- R^2 calculation -----
SS_res1 = sum((y1 - y_fit1).^2);        % Residual sum of squares
SS_tot1 = sum((y1 - mean(y1)).^2);      % Total sum of squares
R21 = 1 - SS_res1/SS_tot1;              % Coefficient of determination

SS_res2 = sum((y2 - y_fit2).^2);        % Residual sum of squares
SS_tot2 = sum((y2 - mean(y2)).^2);      % Total sum of squares
R22 = 1 - SS_res2/SS_tot2;              % Coefficient of determination

SS_res3 = sum((y3 - y_fit3).^2);        % Residual sum of squares
SS_tot3 = sum((y3 - mean(y3)).^2);      % Total sum of squares
R23 = 1 - SS_res3/SS_tot3;              % Coefficient of determination

SS_res4 = sum((y4 - y_fit4).^2);        % Residual sum of squares
SS_tot4 = sum((y4 - mean(y4)).^2);      % Total sum of squares
R24 = 1 - SS_res4/SS_tot4;              % Coefficient of determination
% ---------------------------


% Plot
figure(1);
plot(x, y1, 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(x, y_fit1, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('y1');
legend('Data', 'Least-squares fit');
grid on;

% Display results
fprintf('Slope (m) = %.4f\n', m1);
fprintf('Intercept (c) = %.4f\n', c1);
fprintf('R^2 = %.4f\n', R21);


figure(2);
plot(x, y2, 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(x, y_fit2, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('y2');
legend('Data', 'Least-squares fit');
grid on;

% Display results
fprintf('Slope (m) = %.4f\n', m2);
fprintf('Intercept (c) = %.4f\n', c2);
fprintf('R^2 = %.4f\n', R22);



figure(3);
plot(x, y3, 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(x, y_fit3, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('y3');
legend('Data', 'Least-squares fit');
grid on;

% Display results
fprintf('Slope (m) = %.4f\n', m3);
fprintf('Intercept (c) = %.4f\n', c3);
fprintf('R^2 = %.4f\n', R23);




figure(4);
plot(x, y4, 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(x, y_fit4, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('y4');
legend('Data', 'Least-squares fit');
grid on;

% Display results
fprintf('Slope (m) = %.4f\n', m4);
fprintf('Intercept (c) = %.4f\n', c4);
fprintf('R^2 = %.4f\n', R24);
