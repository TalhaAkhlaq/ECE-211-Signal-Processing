% Talha Akhlaq (Problem Set 6, Problem 3)
clc;
clear;
close all;

% Define constants
R = 1e3;         
C = 10e-9;       
N = 1e4;                    % Number of frequency points for the plot
f = linspace(0, 1e6, N);    % Frequency range from DC to 1 MHz (linear scale)
w = 2*pi*f;     

% Stability constraints and K values
K0 = 4 - 2*sqrt(2);   
Kmax = 4;           
K1 = K0;                 % K1 = K0 (1.1716)
K2 = 0.5*K0 + 0.5*Kmax;  % K2 = 0.5*K0 + 0.5*Kmax (2.5858)
K3 = 0.2*K0 + 0.8*Kmax;  % K3 = 0.2*K0 + 0.8*Kmax (3.4343)

% Natural frequency 
omega_n = sqrt(2 / (R * C));

% PART A: Frequency Response Plots (Linear Scale)
disp('PART A: Frequency Response Plots (Linear Scale)');

% Compute frequency responses 
num_K1 = [K1/(R*C), 0];
den_K1 = [1, (4 - K1)/(R*C), 2/(R^2*C^2)];
[H1, ~] = freqs(num_K1, den_K1, w);

num_K2 = [K2/(R*C), 0];
den_K2 = [1, (4 - K2)/(R*C), 2/(R^2*C^2)];
[H2, ~] = freqs(num_K2, den_K2, w);

num_K3 = [K3/(R*C), 0];
den_K3 = [1, (4 - K3)/(R*C), 2/(R^2*C^2)];
[H3, ~] = freqs(num_K3, den_K3, w);

% Plot magnitude and phase responses 
figure;
subplot(2,1,1);
plot(f/1e3, 20*log10(abs(H1)), 'b', f/1e3, 20*log10(abs(H2)), 'r', f/1e3, 20*log10(abs(H3)), 'g');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Magnitude Response');
legend(['K = ', num2str(K1, '%.4f')], ['K = ', num2str(K2, '%.4f')], ['K = ', num2str(K3, '%.4f')]);
grid on;

subplot(2,1,2);
plot(f/1e3, unwrap(angle(H1))*180/pi, 'b', f/1e3, unwrap(angle(H2))*180/pi, 'r', f/1e3, unwrap(angle(H3))*180/pi, 'g');
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
title('Phase Response');
legend(['K = ', num2str(K1, '%.4f')], ['K = ', num2str(K2, '%.4f')], ['K = ', num2str(K3, '%.4f')]);
grid on;

% PART B: Bode Plots (Logarithmic Frequency Scale)
disp('PART B: Bode Plots (Logarithmic Frequency Scale)');
f_log = logspace(3,6,N);   % Logarithmic frequency axis from 1 kHz to 1 MHz
w_log = 2*pi*f_log;        % Angular frequency for logarithmic scale

% Calculate frequency responses
[H1_log, ~] = freqs(num_K1, den_K1, w_log);
[H2_log, ~] = freqs(num_K2, den_K2, w_log);
[H3_log, ~] = freqs(num_K3, den_K3, w_log);

% Plot Bode magnitude and phase responses
figure;
subplot(2,1,1);
semilogx(f_log/1e3, 20*log10(abs(H1_log)), 'b', f_log/1e3, 20*log10(abs(H2_log)), 'r', f_log/1e3, 20*log10(abs(H3_log)), 'g');
xlabel('Frequency (kHz)');
ylabel('Magnitude (dB)');
title('Bode Magnitude Response');
legend(['K = ', num2str(K1, '%.4f')], ['K = ', num2str(K2, '%.4f')], ['K = ', num2str(K3, '%.4f')]);
grid on;

subplot(2,1,2);
semilogx(f_log/1e3, unwrap(angle(H1_log))*180/pi, 'b', f_log/1e3, unwrap(angle(H2_log))*180/pi, 'r', f_log/1e3, unwrap(angle(H3_log))*180/pi, 'g');
xlabel('Frequency (kHz)');
ylabel('Phase (degrees)');
title('Bode Phase Response');
legend(['K = ', num2str(K1, '%.4f')], ['K = ', num2str(K2, '%.4f')], ['K = ', num2str(K3, '%.4f')]);
grid on;
