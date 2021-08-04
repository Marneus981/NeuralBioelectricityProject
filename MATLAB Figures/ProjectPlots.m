%% Transmembrane Voltage Vm
% Time
Vr = -60e-3;
A = 45e-3; Su = 2.2; tu = 4;
B = 45e-3; Sd = .73; td = 8;

Vm_t = @(t) [Vr + A.*tanh(Su.*(t-tu)) - B.*tanh(Sd.*(t-td))];

% Space (for RIGHT moving wave i.e. (-))
u = 2; 
x0 = 7;
t0 = 1.7;

% Transmembrane Voltage and Derivatives in T and X
Vm_t_x = @(t,x) [Vr + A.*tanh(Su.*(t-tu - 1/u.*(x-x0))) - B.*tanh(Sd.*(t-td - 1/u.*(x-x0)))];
dVmdx_t_x = @(t,x) [-49.5e-3 .* (sech(Su.*(t - tu - 1/u.*(x-x0)))).^2 + 16.425e-3 .*(sech(Sd*(t - td - 1/u.*(x-x0)))).^2]; 
ddVmddx_t_x = @(t,x) [-108.9e-3 .* (sech(Su.*(t - tu - 1/u.*(x-x0)))).^2 .* tanh(Su*(t - tu - 1/u.*(x-x0))) + 11.99e-3 .* (sech(Sd.*(t - td - 1/u.*(x-x0))).^2).*tanh(Sd.*(t - td - 1/u.*(x-x0)))];

% Time Plot
T = 0:0.01:20;
Vm_T = Vm_t(T).*1000;

figure('name', 'Transmembrane Voltage - Time');
plot(T, Vm_T);
grid MINOR;
title("\textbf{Transmembrane Voltage - Time at} $x_{0}$ = 7 mm", 'interpreter', 'latex');
ylabel('$V_{m}$[mV]', 'interpreter', 'latex');
xlabel('\textbf{Time} [ms]', 'interpreter', 'latex');
%legend('$x_{0}$ = 7mm','interpreter','latex');

% Space Plot
X = -20:0.01:20;
Vm_X = Vm_t_x(t0, X).*1000;

figure('name', 'Transmembrane Voltage - Space');
plot(X, Vm_X);
grid MINOR;
title("\textbf{Transmembrane Voltage - Space at} $t_{0}$ = 1.7 ms", 'interpreter', 'latex');
ylabel('$V_{m}$[mV]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');
%legend('$t_{0}$ = 1.7ms','interpreter','latex');

% 3D Space-Time Plot
[Xmesh, Tmesh] = meshgrid(-20:0.5:20,0:0.5:20);
figure('name', 'Transmembrane Voltage - SpaceTime');
surf(Xmesh, Tmesh, 1000.*Vm_t_x(Tmesh,Xmesh));
grid MINOR;
title("\textbf{Transmembrane Voltage - SpaceTime}", 'interpreter', 'latex');
zlabel('$V_{m}$[mV]', 'interpreter', 'latex');
ylabel('\textbf{Time} [ms]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');

%% Transmembrane Current Im
% Givens
a = 100e-6; % micrometers
sigma_i = 0.02/1e-2; % Siemens per centimeters
sigma_e = inf; % Siemens per centimeters

im = @(x) [pi*a.^2*sigma_i*ddVmddx_t_x(t0, x)];
im_t_x = @(t,x) [pi*a.^2*sigma_i*ddVmddx_t_x(t, x)];

% Time Plot
T = 0:0.01:20;
im_T_x0 = im_t_x(T,x0)*10e9;

figure('name', 'Transmembrane Current - Time');
plot(T, im_T_x0);
grid MINOR;
title("\textbf{Transmembrane Current - Time at} $x_{0}$ = 7 mm", 'interpreter', 'latex');
ylabel('$\i_{m}$ [nA]', 'interpreter', 'latex');
xlabel('\textbf{Time} [ms]', 'interpreter', 'latex');
%legend('$x_{0}$ = 7m','interpreter','latex');

% Space Plot
X = -20:0.01:20;
im_t0_X = im_t_x(t0,X)*10e9;

figure('name', 'Transmembrane Current - Space');
plot(X, im_t0_X);
grid MINOR;
title("\textbf{Transmembrane Current - Space at} $t_{0}$ = 1.7 ms", 'interpreter', 'latex');
ylabel('$\i_{m}$ [nA]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');
%legend('$t_{0}$ = 1.7s','interpreter','latex');

% 3D Space-Time Plot
[Xmesh, Tmesh] = meshgrid(-20:0.5:20,0:0.5:20);
figure('name', 'Transmembrane Current - SpaceTime');
surf(Xmesh, Tmesh, 10e9.*im_t_x(Tmesh,Xmesh));
grid MINOR;
title("\textbf{Transmembrane Current - SpaceTime}", 'interpreter', 'latex');
zlabel('$\i_{m}$ [nA]', 'interpreter', 'latex');
ylabel('\textbf{Time} [ms]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');

%% Lumped Monopoles
% Givens
% x from -20 to 20
x1 = fzero(im,-20);
x2 = fzero(im, 0);
x3 = fzero(im, 20);

% Monopole Current
M1 = pi*a^2*sigma_i*(dVmdx_t_x(t0, x1) - dVmdx_t_x(t0, -20));
M2 = pi*a^2*sigma_i*(dVmdx_t_x(t0, x2) - dVmdx_t_x(t0, x1));
M3 = pi*a^2*sigma_i*(dVmdx_t_x(t0, x3) - dVmdx_t_x(t0, x2));
M4 = pi*a^2*sigma_i*(dVmdx_t_x(t0, 20) - dVmdx_t_x(t0, x3));

% Monopole Position
imx = @(x) im(x).*x;
o1 = 1/M1*integral(imx, -20, x1);
o2 = 1/M2*integral(imx, x1, x2);
o3 = 1/M3*integral(imx, x2, x3);
o4 = 1/M4*integral(imx, x3, 20);

%% Extracellular Potentials 
y = 2; 
c = 7.958e-3;

% Space
r = @(o, x) (sqrt((x-o).^2+1));
phi_e_x = @(x) [c.*(M1./r(o1,x) + M2./r(o2,x) + M3./r(o3,x) + M4./r(o4,x))];

% Time
phi_e_t = @(t) [phi_e_x(x0 - u.*(t-t0))];

% Space-Time
phi_e_t_x = @(t,x) [phi_e_x(x - u.*(t-t0))];

% Time Plot
T = 0:0.01:20;
phi_e_T = 10e12.*phi_e_t(T);

figure('name', 'Extracellular Potential - Time');
plot(T, phi_e_T);
grid MINOR;
title("\textbf{Extracellular Potential - Time at} $x_{0}$ = 7mm", 'interpreter', 'latex');
ylabel('$\Phi_{e}$ [pV]', 'interpreter', 'latex');
xlabel('\textbf{Time} [ms]', 'interpreter', 'latex');
%legend('$x_{0}$ = 7m','interpreter','latex');

% Space Plot
X = -20:0.01:20;
phi_e_X = 10e12.*phi_e_x(X);

figure('name', 'Extracellular Potential - Space');
plot(X, phi_e_X);
grid MINOR;
title("\textbf{Extracellular Potential - Space at} $t_{0}$ = 1.7ms", 'interpreter', 'latex');
ylabel('$\Phi_{e}$ [pV]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');
%legend('$t_{0}$ = 1.7s','interpreter','latex');

% 3D Space-Time Plot
[Xmesh, Tmesh] = meshgrid(-20:0.5:20,0:0.5:20);
figure('name', 'Extracellular Potential - SpaceTime');
surf(Xmesh, Tmesh, 10e12.*phi_e_t_x(Tmesh,Xmesh));
grid MINOR;
title("\textbf{Extracellular Potential - SpaceTime}", 'interpreter', 'latex');
zlabel('$\Phi_{e}$ [pV]', 'interpreter', 'latex');
ylabel('\textbf{Time} [ms]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');

%% Lumped Monopole I_m and Phi_e plots
% Extracellular Potential with Monopoles
figure('name', 'Extracellular Potential with Monopoles - Space');
hold on;

scatter(o1, 0, 5*10e9*abs(M1), 'b','filled');
scatter(o2, 0, 5*10e9*abs(M2), 'r','filled');
scatter(o3, 0, 5*10e9*abs(M3), 'b','filled'); 
q1 = quiver(o1,0,0, 10e12.*phi_e_x(o1), 1, 'ShowArrowHead', 'off', 'color', 'b');
q2 = quiver(o2,0,0, 10e12.*phi_e_x(o2), 1, 'ShowArrowHead', 'off', 'color', 'r');
q3 = quiver(o3,0,0, 10e12.*phi_e_x(o3), 1, 'ShowArrowHead', 'off', 'color', 'b');


plot(X, phi_e_X, 'k');

grid MINOR;
title("\textbf{Extracellular Potential with Monopoles - Space at} $t_{0}$ = 1.7 ms", 'interpreter', 'latex');
ylabel('$\Phi_{e}$ [pV]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');
legend('$M_{1}$ = ' + string(M1*10e9) + " nA", '$M_{2}$ = ' + string(M2*10e9) + " nA", '$M_{3}$ = ' + string(M3*10e9) + " nA",'interpreter','latex');
hold off; 

% Transmembrane Current with Monopoles
figure('name', 'Transmembrane Current with Monopoles - Space');
hold on;
aplot1 = area(X(1:1441), im_t0_X(1:1441), 'FaceColor', 'b');
aplot2 = area(X(1441:2241), im_t0_X(1441:2241), 'FaceColor', 'r');
aplot3 = area(X(2241:end), im_t0_X(2241:end), 'FaceColor', 'b');

bar(o1, 10e9*M1, 0.01, 'k');
scatter(o1, 10e9*M1, 20, 'b');

bar(o2, 10e9*M2, 0.01, 'k');
scatter(o2, 10e9*M2, 20, 'r');


bar(o3, 10e9*M3, 0.01, 'k');
scatter(o3, 10e9*M3, 20, 'b'); 

%bar([o1, o2, o3], [10e9*M1, 10e9*M2, 10e9*M3], 0.01, 'b');
%scatter([o1, o2, o3], [10e9*M1, 10e9*M2, 10e9*M3], [20 20 20], [1 0 0; 1 0 0; 1 0 0]);
grid MINOR;
title("\textbf{Transmembrane Current with Monopoles - Space at} $t_{0}$ = 1.7 ms", 'interpreter', 'latex');
ylabel('$\i_{m}$ [nA]', 'interpreter', 'latex');
xlabel('\textbf{x} [mm]', 'interpreter', 'latex');
legend('Sources','Sinks','interpreter','latex');
hold off; 

