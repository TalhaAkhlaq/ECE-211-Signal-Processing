%% Talha Akhlaq (Problem Set 8)
clc; clear; close all;

%% Part 1: Heavy‐Tail Distributions

% 1(a) Generate N = 10^6 iid samples
N        = 1e6;                                 % total number of samples
x_gauss  = randn(N,1);                          % Gaussian N(0,1)
x_t      = trnd(5, N,1) * sqrt((5-2)/5);        % Student’s t (ν=5), scaled to Var=1
x_cauchy = 0.544 * tan(pi*(rand(N,1)-0.5));     % Cauchy(scale=α=0.544)

% 1(b) Fundamental assumption
%     Law of Large Numbers ⇒ empirical frequency → true probability

% 1(c) Estimate tail probabilities P(|X|>2)
p_gauss  = mean(abs(x_gauss)  > 2);             % P(|Gaussian|>2)
p_t      = mean(abs(x_t)      > 2);             % P(|Student’s t|>2)
p_cauchy = mean(abs(x_cauchy) > 2);             % P(|Cauchy|>2)

% 1(c) Plot each sequence with ±2 thresholds
figure(1);
plot(x_gauss, '.', 'MarkerSize', 1); hold on;   % raw Gaussian
yline( 2, 'r--'); yline(-2, 'r--');              
grid on; title('Gaussian samples');        
xlabel('n'); ylabel('x[n]');

figure(2);
plot(x_t, '.', 'MarkerSize', 1); hold on;
yline(2, 'r--'); yline(-2, 'r--');
grid on;
title('Student''s t samples');
xlabel('n'); ylabel('x[n]');

figure(3);
plot(x_cauchy, '.', 'MarkerSize', 1); hold on;
yline(2, 'r--'); yline(-2, 'r--');
grid on;
title('Cauchy (\alpha = 0.544, as per PDF formula) Data');
xlabel('Sample Index');
ylabel('Value');

% 1(c) Reorganize into 10 segments of 100,000 each and compute means
segLen   = 1e5;                                 % segment length
nSeg     = N / segLen;                         % number of segments
mu_g     = zeros(nSeg,1);
mu_t2    = zeros(nSeg,1);
mu_c2    = zeros(nSeg,1);
for k = 1:nSeg
    idx       = (k-1)*segLen + (1:segLen);
    mu_g(k)   = mean(x_gauss(idx));             % Gaussian block mean
    mu_t2(k)  = mean(x_t(idx));                 % Student’s t block mean
    mu_c2(k)  = mean(x_cauchy(idx));            % Cauchy block mean
end

figure(4);

subplot(3,1,1);
stem(1:nSeg, mu_g, 'filled'); grid on;
xlabel('Block Index'); ylabel('Mean');
title('Block Means – Gaussian');

subplot(3,1,2);
stem(1:nSeg, mu_t2, 'filled'); grid on;
xlabel('Block Index'); ylabel('Mean');
title('Block Means – Student''s t (ν=5)');

subplot(3,1,3);
stem(1:nSeg, mu_c2, 'filled'); grid on;
xlabel('Block Index'); ylabel('Mean');
title('Block Means – Cauchy (α=0.544)');

% Observation (1c): Cauchy block means vary widely; sample mean not reliable

%% Part 2: ARMA and AR Modeling

% Model: x[n] = v[n] + 0.4 v[n-1] + 0.2 v[n-2] + 1.6 x[n-1] − 0.81 x[n-2]
b = [1, 0.4, 0.2];    % MA(2) coefficients
a = [1, -1.6, 0.81];  % AR(2) coefficients

% 2(a) Pole–zero analysis (minimum‐phase check)
[z,p,~] = tf2zp(b,a);                         % compute zeros and poles
theta   = linspace(0,2*pi,300);
figure(5);
plot(cos(theta), sin(theta), 'w--'); hold on; % white unit circle
plot(real(z), imag(z), 'bo', 'MarkerSize',8); % zeros in blue
plot(real(p), imag(p), 'rx', 'MarkerSize',8); % poles in red
axis equal; grid on; grid minor;
set(gca, ...                                   % black background
    'Color','k','XColor','w','YColor','w', ...
    'GridColor','w','MinorGridColor','w');
title('Pole–Zero Plot of H(z)','Color','w');
xlabel('Real','Color','w'); ylabel('Imag','Color','w');
% Confirm: all poles & zeros inside unit circle → minimum‐phase

% 2(b) Generate ARMA process & estimate autocorrelation r_x[m]
N2    = 1e5; sigma2 = 2;
v2    = sqrt(sigma2) * randn(N2,1);             % white‐noise input
x2    = filter(b, a, v2);                       % ARMA output
M     = 6;                                      % max lag
r_x   = zeros(M+1,1);
for m = 0:M
    r_x(m+1) = x2(1:N2-m)' * x2(1+m:N2) / (N2-m);% unbiased estimate
end
r_full = [r_x(end:-1:2); r_x]; lags = -M:M;

figure(6);
stem(lags, r_full, 'filled'); grid on;
title('Autocorrelation r_x[m]');
xlabel('Lag m'); ylabel('r_x[m]');

% 2(c) PSD estimation via Welch’s method; mark peak ω₀
[Px,w]  = pwelch(x2, hamming(512), 256, 512);
[~,i0]  = max(Px); omega0 = w(i0);               % peak frequency

figure(7);
plot(w, Px, 'LineWidth',1.5); hold on;
xline(omega0, 'r--'); grid on;
title('Welch PSD Estimate');
xlabel('\omega (rad/sample)'); ylabel('S_x(\omega)');
legend('PSD','\omega_0','Location','best');
% ω₀ corresponds to spectral peak due to filter poles

% 2(d) AR(4) via Yule–Walker; final comparison stem plot
[a4, var4] = aryule(x2, 4);                     % AR(4) coefficients + var
x0         = filter(1, a4, v2);                 % AR(4) reconstruction

figure(8);
h1 = stem(1:100, x2(1:100), 'b-o','LineWidth',1.2); hold on;
h2 = stem(1:100, x0(1:100), 'r--x','LineWidth',1.2);
grid on;
title('Comparison of First 100 Samples','FontWeight','normal');
xlabel('Sample Index n'); ylabel('Amplitude');
legend([h1 h2], ...
    'x[n] (Original ARMA(2,2))', ...
    'x_0[n] (AR(4) model output)', ...
    'Location','best');

% Estimate autocorrelation of x0 (AR(4) output)
r0 = zeros(2*M+1,1);
for m = -M:M
    if m < 0
        r0(M+1+m) = x0(1:end+m)' * x0(1-m:end) / (N2 - abs(m));
    else
        r0(M+1+m) = x0(1+m:end)' * x0(1:end-m) / (N2 - abs(m));
    end
end

% Plot original and AR(4) autocorrelations together
figure(9);
stem(-M:M, r_full, 'bo', 'filled'); hold on;
stem(-M:M, r0, 'r--o');
grid on;
xlabel('Lag');
ylabel('Correlation Value');
title('Comparison of Autocorrelations');
legend('ARMA(2,2)', 'AR(4) Model Output', 'Location', 'Best');

% Confirm: x[n] and x_0[n] differ pointwise but share stats

% 2(b.4) Build Toeplitz matrix R1 from r_x
R1 = toeplitz(r_x(1:7));
disp('R1 ='); disp(R1);

% 2(b.5) Compute eigenvalues
eigenvals = eig(R1);
disp('Eigenvalues of R1:'); disp(eigenvals);

% 2(b.6) Estimate R2 from raw data
X = zeros(N2 - M, M + 1);
for i = 1:(M + 1)
    X(:, i) = x2(i:(end - M + i - 1));
end
R2 = (X' * X) / size(X,1);

% Compute norm difference
R_diff_norm = norm(R1 - R2);
disp('‖R1 − R2‖ ='); disp(R_diff_norm);