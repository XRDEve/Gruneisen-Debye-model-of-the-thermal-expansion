% DIFFRACTION DATA ANALYSIS % DIFFRACTION DATA ANALYSIS 
% CODE INFORMATION 
% Fits the 2nd order Gruneissen approximation model to experimental volume data and calculates the
% volumetric thermal expansion aV(T), etc
% Version 3.0 2018
% @E.K.Tripoliti ucfbetr@ucl.ac.uk
%
rng(9845, 'twister')                                                  
%
% INPUT VARIABLES
kB = 1.38064852e-23;            % Bolzmann's constant J/K
N = input('Atoms in the unit cell:');                                 
gth = input('Gruneisen Parameter (if not known set equal to 1):');       
n = input('Formula Units:');                                                                                  
AvNum = 6.0221409e+23;         % Avogadro's Number
T_exp =data(1:end,1);                                                       
V_exp =data(1:end,2);  SV = data(1:end,3);                                 

%% INTEGRAL OF INTERNAL ENERGY 
% Using the K24 approximation of Khishchenko (2020) %% 2nd ORDER GRUNEISEN-DEBYE MODEL SOLUTION
A03 = (pi^4)/5;
A10 = ((8*pi^4)/15) - 49;
A11 = ((8*pi^4)/5) - (219/2) - 12*A10;
A12 = ((16*pi^4)/5) - 117 - 36*A10 -8*A11;
A13 = ((16*pi^4)/5) - 39 - 24*A10 -12*A11 - 4*A12;
A20 = ((4*pi^4)/15) - 1 - A10 - A11 - (A12/2) - A13/6;
A21 = ((2*pi^4)/5) - A11 - A12 - A13/2;
A22 = ((2*pi^4)/5) - A12 - A13;
A23 = ((pi^4)/5) - A13;
%% Set model 
% Initial model parameters
% x = [V0 ThD Q b]; Initial model 'x' and bounds 'xl' and 'xu' can be
% changed if they need to be more constrained
x = [288 200 1 0];
xl =  [280  600 1 1];
xu = [1000  1000 9 10];
% Optimization 
options = optimoptions( 'lsqnonlin','MaxIterations', 1000, 'FunctionTolerance', 1e-31, 'Algorithm',...
'levenberg-marquardt', 'Display', 'final-detailed', 'FiniteDifferenceType', 'central', 'MaxFunctionEvaluations', ...
4000, 'OptimalityTolerance', 1e-31, 'StepTolerance', 1e-31 );
% Know that the levenberg-marquardt algorithm requires bounds. If bounds
% are not needed then use trust-region-reflective.
% Creat anonymous function 
fun = @(x) (((x(1) + ((x(1).*(9*kB*N.*T_exp.*((T_exp./x(2)).^3).*(A03./((x(2)./T_exp).^(3)) - ...
(A10 + A11./(x(2)./T_exp) + A12./((x(2)./T_exp).^2) + A13./((x(2)./T_exp).^3)) .* exp(-1.*(x(2)./T_exp)) - ...
(A20 + A21./(x(2)./T_exp)+ A22./((x(2)./T_exp).^2) +...
A23./((x(2)./T_exp).^3)) .* exp(-2.*(x(2)./T_exp)))./(3./((x(2)./T_exp).^3))))./((x(3)*1e-17)-x(4).*(9*kB*N.*T_exp.*((T_exp./x(2)).^3).*(A03./((x(2)./T_exp).^(3)) ...
- (A10 + A11./(x(2)./T_exp) + A12./((x(2)./T_exp).^2) + A13./((x(2)./T_exp).^3)) .* exp(-1.*(x(2)./T_exp)) -  (A20 + A21./(x(2)./T_exp) + A22./((x(2)./T_exp).^2) ...
+ A23./((x(2)./T_exp).^3)) .* exp(-2.*(x(2)./T_exp)))./(3./((x(2)./T_exp).^3)))))-V_exp).^2));
% Non-linear least squares algorithm
[x,resnorm] = lsqnonlin(fun,x,xl,xu,options);
% Export fitted values 
display([' V0:'     num2str(x(1))]);                display([' TH_D:'   num2str(x(2))]);
display([' Q:'      num2str(x(3))]);                 display([' b:'      num2str(x(4))]);
%% CALCULATIONS 
% Solve internal energy with correct Debye Temperature 
DT = (x(2)./T_exp);                                                                                                                                                                           % ThD/T
K24 = (A03./(DT.^(3)) - (A10 + A11./DT + A12./(DT.^2) + A13./(DT.^3)) .* exp(-1.*DT) - ...                                                          % K24 approximation of U(T) integral 
(A20 + A21./DT + A22./(DT.^2) + A23./(DT.^3)) .* exp(-2.*DT))./(3./(DT.^3));
U_T = (9*kB*N.*T_exp.*((T_exp./x(2)).^3).*K24);                                                                                                                              % Internal energy U(T)
% Unit-cell Volume 
V_T = x(1) + ((x(1).* (9*kB*N.*T_exp.*((T_exp./x(2)).^3).*K24))./((x(3)*1e-17)-(x(4).* (9*kB*N.*T_exp.*((T_exp./x(2)).^3).*K24))));   % Calculated V(T) with V0,thD,Q and b values calculated from the model 
% Quality of fit
SR = (V_exp - V_T).^2; SSR = sum(SR);                                                                                                                                           % Sum of squared errors
display([' Sum of Standard Differences:'     num2str(SSR)]);  

%% EXTRAPOLATIONS 
T_pred = (4:0.2:1200);  % T_pred can change. Extreme extrapolations are not adviced.                                                                                                                                                                    % Extrapolated Temperature; this can change
T_pred = T_pred(:);                                                                                                                                                                          % Convert T_exp into a n:1 array
DT_pred = (x(2)./T_pred);                                                                                                                                                                % Extrapolated ThD/T
K24_pred = (A03./(DT_pred.^(3)) - (A10 + A11./DT_pred + A12./(DT_pred.^2) + A13./(DT_pred.^3)) .* exp(-1.*DT_pred) - ...    % Extrapolated K24
(A20 + A21./DT_pred + A22./(DT_pred.^2) + A23./(DT_pred.^3)) .* exp(-2.*DT_pred))./(3./(DT_pred.^3));
U_T_pred = (9*kB*N.*T_pred.*((T_pred./x(2)).^3).*K24_pred);                                                                                                        % Extrapolated U(T)               
V_T_pred = x(1) + ((x(1).* (9*kB*N.*T_pred.*((T_pred./x(2)).^3).*K24_pred))./((x(3)*1e-17)-...                                                        % Extrapolated V(T)
(x(4).* (9*kB*N.*T_pred.*((T_pred./x(2)).^3).*K24_pred))));

%% PLOTTING
%No.1: Unit-cell Volume vs Temperature & differences between experimental data and model
figure(1)
% plot experimenta values
subplot(2,2,1)
plot(T_exp, V_exp,  'ko','MarkerSize',5,'LineWidth',0.8,'MarkerFaceColor', 'w')
hold on; 
xlim([0 max(T_pred)])
% plot extrapolated model
plot(T_pred, V_T_pred, 'k-','LineWidth',0.5)
xlabel('T(K)', 'FontWeight','bold', 'fontsize', 12)
ylabel('V(T) (Å^3)', 'FontWeight','bold', 'fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off
% plot differences
subplot(2,2,2)
errorbar(T_exp, (V_exp - V_T), SV,  'ko','MarkerSize',5,'LineWidth',0.8,'MarkerFaceColor', 'w')
hold on; 
xlim([0 max(T_pred)])
linex=[0 max(T_pred)];
liney = [0 0];
plot(linex, liney,  ':','LineWidth',1,'Color', [0.8 0.8 0.8])
xlabel('T (K)', 'FontWeight','bold', 'fontsize', 12)
ylabel('V(T)_o_b_s - V(T)_c_a_l_c (Å^3)', 'FontWeight','bold', 'fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off
% 
%% STANDARD ERRORS
% Partial Derivatives and Weights 
D1 = ((x(1)*(1+0.001).*U_T./((x(3).*(1e-17))-x(4).*U_T) + x(1)*(1+0.001))-V_T)./(x(1)*0.001);
DoverT = x(2).*(1+0.001)./T_exp;
K24_1 = (A03./(DoverT.^(3)) - (A10 + A11./DoverT + A12./(DoverT.^2) + A13./(DoverT.^3)) .* exp(-1.*DoverT) - ...
    (A20 + A21./DoverT + A22./(DoverT.^2) + A23./(DoverT.^3)) .* exp(-2.*DoverT))./(3./(DoverT.^3));
UT_1 = (9*kB*N.*T_exp.*((T_exp./(x(2).*(1+0.001))).^3).*K24_1);
D2 = (((((x(1).*UT_1 )./((x(3)*1e-17)-x(4).*UT_1))) +x(1)) - V_T)/(x(2)*0.001);
D3 = (((((x(1).*UT_1)./(((x(3)*1e-17)*(1+0.001))-(x(4).*UT_1)))) +x(1)) - V_T)/ ((x(3)*1e-17)*0.001);
D4 = ((((x(1).*UT_1)./((x(3)*1e-17)-(x(4)*(1+0.001)).*UT_1))+x(1))-V_T)./(x(4)*0.001);
% For weighted error use SV = experimental, for equally weighted use SV = 1.
D11 = SV.*(D1.^2);           SD11 = sum(D11,'all');
D12 = SV.*D1.*D2;            SD12 = sum(D12,'all');
D13 = SV.*D1.*D3;            SD13 = sum(D13,'all');
D14 = SV.*D1.*D4;            SD14 = sum(D14,'all');
%
D22 = SV.*(D2.^2);           SD22 = sum(D22,'all');
D23 = SV.*D2.*D3;            SD23 = sum(D23,'all');
D24 = SV.*D2.*D4;            SD24 = sum(D24,'all');
%
D33 = SV.*(D3.^2);           SD33 = sum(D33,'all');
D34 = SV.*D3.*D4;            SD34 = sum(D34,'all');
%
D44 = SV.*(D4.^2);           SD44 = sum(D44,'all');
% aij & bij matrices
aij = zeros(4,4);
aij(1,1) = SD11; aij(1,2) = SD12; aij(1,3) = SD13; aij(1,4) = SD14;
aij(2,1) = SD12; aij(2,2) = SD22; aij(2,3) = SD23; aij(2,4) = SD24; 
aij(3,1) = SD13; aij(3,2) = SD23; aij(3,3) = SD33; aij(3,4) = SD34;
aij(4,1) = SD14; aij(4,2) = SD24; aij(4,3) = SD34; aij(4,4) = SD44;
bij = inv(aij);                         % Inversion matrix
Cij = corrcoef(bij);                    % Correlation matrix 
NoU = 4;                                % Number of unknown variables in the model since it is the 2nd-order
% e.s.d 
Error_V0 = sqrt(bij(1,1).*SSR./(length(data) - NoU));
Error_ThD = sqrt(bij(2,2).*SSR./(length(data) - NoU));
Error_Q = sqrt(bij(3,3).*SSR./(length(data) - NoU));
Error_b = sqrt(bij(4,4).*SSR./(length(data) - NoU));
% display errors in the command window
display([' Error V0:'     num2str(Error_V0)]);             display([' Error TH_D:'   num2str(Error_ThD)]);
display([' Error Q:'      num2str(Error_Q)]);              display([' Error b:'      num2str(Error_b)]);
%% THERMAL EXPANSION
% a_V(T) =1/T*[diff(V)/diff(T)] 
alpha_calc = (1./V_T(2:end)).*(diff(V_T)./diff(T_exp));               % aV(T) from experimenta data. This is a point-by-point differentiation
alpha_exp = (1./V_exp(2:end)).*(diff(V_exp)./diff(T_exp));            % aV(T) from model 
alpha_pred = (1./V_T_pred(2:end)).*(diff(V_T_pred)./diff(T_pred));    % aV(T) from extrapolated model between 4 and 1500 K; T_pred can change in the previous paragraph 
% No.2: Thermal Expansion vs Temperature (observed with model fit)
figure(2)
% plot aV(T) from experimental data
subplot(2,2,1)
plot(T_exp(2:end), alpha_exp,  'ko',  'Markersize', 5,'LineWidth',0.8)
hold on;
% plot aV(T) from extrapolated model as line fitted to the data
plot(T_pred(2:end), alpha_pred,  'k-','LineWidth',0.5)
xlim([0 max(T_pred)])
ylim([0 max(alpha_exp)])
xlabel('T(K)', 'FontWeight','bold', 'fontsize', 12)
ylabel('\alpha_V(T) (K^-^1)', 'FontWeight','bold', 'fontsize', 12)
set(gca,'FontSize', 13, 'FontWeight','bold')
hold off
% plot dofferences
subplot(2,2,2)
plot(T_exp(2:end), (alpha_exp-alpha_calc), 'ko',  'Markersize', 5,'LineWidth',0.8,'MarkerFaceColor', 'w')
hold on;
ylim([-1e-05 1e-05])
xlim([0 max(T_pred)])
linex=[0 max(T_pred)];
liney = [0 0];
plot(linex, liney,  ':','LineWidth',1,'Color', [0.8 0.8 0.8])
xlabel('T (K)', 'FontWeight','bold', 'fontsize', 12)
ylabel(' \alpha_V(T)_o_b_s - \alpha_V(T)_c_a_l_c (K^-^1)', 'FontWeight','bold','fontsize', 12)
set(gca,'FontSize', 13, 'FontWeight','bold')
hold off
%% ADDITIONAL THERMAL QUANTITIES
% No.3: The Internal Energy vs Temperature
figure(3)
plot(T_exp, U_T, 'ko',  'Markersize', 5,'LineWidth',0.8,'MarkerFaceColor', 'w')
hold on
plot(T_pred, U_T_pred,  'k-','LineWidth',0.5)
xlabel('T (K)', 'FontWeight','bold', 'fontsize', 12)
ylabel('U(T) (J/K)', 'FontWeight','bold', 'fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
xlim([0 max(T_pred)])
hold off
% Heat Capacities 
% Solving the Integral 
% Model
z = x(2)./T_pred;
fun1 =@(z)((z.^4).*exp(z))./((exp(z)-1).^2);
fmodel=[];
for i=1:length(z)
fmodel(i) = integral(fun1, 0, z(i));
i=i+1;
end
fmodel = fmodel(:);
% Experiment
z1 = x(2)./T_exp;
fun2 =@(z1)((z1.^4).*exp(z1))./((exp(z1)-1).^2);
fexperiment=[];
for i=1:length(z1)
fexperiment(i) = integral(fun2, 0, z1(i));
i=i+1;
end
fexperiment = fexperiment(:);
% Calculating CV for experimental and extrapolated values 
CV_exp = 9*(n)*AvNum*kB.*((T_exp./x(2)).^3).*fexperiment;
CV_pred = 9*(n)*AvNum*kB.*((T_pred./x(2)).^3).*fmodel;
% No.4: Heat Capacity vs Temperature
figure(4)
subplot(3,1,1)
plot(T_exp, CV_exp, 'ko',  'Markersize', 5,'LineWidth',0.8,'MarkerFaceColor', 'w')
hold on
plot(T_pred, CV_pred,  'k-','LineWidth',0.5)
xlabel('T (K)', 'FontWeight','bold', 'fontsize', 12)
ylabel('C_V (J mol^-^1 K^-^1)', 'FontWeight','bold', 'fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off
% 1st Derivative of the Incompressibility K'0
K0_prime = 2*x(4) + 1;
display([' K0_prime :' num2str((K0_prime))]);
% Calculating CP for experimental and extrapolated values: gamma needs to be known. If not it is set equal to 1.
CP_exp = CV_exp(2:end)./(1+(T_exp(2:end).*gth.*alpha_exp));
CP_pred = CV_pred(2:end)./(1+(T_pred(2:end).*gth.*alpha_pred));
% plot CP vs Temperature
subplot(3,1,2)
plot(T_exp(2:end), CP_exp, 'ko',  'Markersize', 5,'LineWidth',0.8,'MarkerFaceColor', 'w')
hold on
plot(T_pred(2:end), CP_pred,  'k-','LineWidth',0.5)
xlabel('T (K)', 'FontWeight','bold', 'fontsize', 12)
ylabel('C_P (J mol^-^1 K^-^1)', 'FontWeight','bold', 'fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off
% plot CV-CP
subplot(3,1,3)
plot(T_exp(2:end), CV_exp(2:end) - CP_exp, 'ko',  'Markersize', 5,'LineWidth',0.8,'MarkerFaceColor', 'w')
hold on
plot(T_pred(2:end), CV_pred(2:end)-CP_pred,  'k-','LineWidth',0.5)
xlabel('T (K)', 'FontWeight','bold', 'fontsize', 12)
ylabel('C_V - C_P (J mol^-^1 K^-^1)', 'FontWeight','bold', 'fontsize', 12)
set(gca,'fontsize', 13,  'FontWeight','bold')
hold off
%








