%% Talha Akhlaq (Problem Set 7)
clc; clear; close all;

%% Part I: Signal Generation
% Define dimensions, sparsity, snapshot count, and power levels
M     = 100;         % ambient dimension of signal space
K     = 20;          % number of nonzero entries per source
L     = 3;           % number of underlying sources
N     = 200;         % number of snapshots (time samples)
PdB   = [0, -2, -4]; % source power levels in dB
Pn_dB = 10;          % noise power in dB

% Generate source vectors S, source coefficients B, noise V, and measurements A
[S, B, V, A] = generateSignal(M, N, K, PdB, Pn_dB);
R = (A * A') / N;    % sample correlation matrix

%% Part II: Analysis
% Compute SVD of A for signal subspace
[U, Sigma, ~] = svd(A);
sing_vals = diag(Sigma);

% Eigen-decomposition of R for comparison
[EigVec, EigValM] = eig(R);
[eig_vals, idx]   = sort(diag(EigValM), 'descend');
eig_vecs = EigVec(:, idx);

% Plot singular values
figure;
stem(sing_vals, 'filled');
title('Singular values σ_k of A  (M=100, N=200)');
xlabel('Singular-value index k');
ylabel('Singular value σ_k');
xline(3, '--', 'L=3', 'LabelOrientation','horizontal');
ylim([0, max(sing_vals)*1.1]);
xticks(0:10:length(sing_vals));
yticks(0:2:ceil(max(sing_vals)*1.1));
grid on; box on;

% Plot eigenvalues
figure;
stem(eig_vals, 'filled');
title('Eigenvalues λ_k of R  (M=100, N=200)');
xlabel('Eigenvalue index k');
ylabel('Eigenvalue λ_k');
ylim([0, max(eig_vals)*1.1]);
xticks(0:10:M);
yticks(0:0.2:ceil(max(eig_vals)*1.1/0.2)*0.2);
grid on; box on;

% Display ratio between 3rd and 4th modes
sigma_ratio  = sing_vals(3)/sing_vals(4);
lambda_ratio = eig_vals(3)/eig_vals(4);
fprintf('σ₃/σ₄ = %.3f\n', sigma_ratio);
fprintf('λ₃/λ₄ = %.3f\n\n', lambda_ratio);

% Build projector onto noise subspace and invert R
UL      = U(:,1:L);
P_noise = eye(M) - UL*UL';
R_inv   = inv(R);

% MUSIC & MVDR at true sources and random directions
numTest = 20;
SM_true = zeros(L,1);
SM_rand = zeros(numTest,1);
MV_true = zeros(L,1);
MV_rand = zeros(numTest,1);

for i = 1:L
    s_i        = S(:,i);
    SM_true(i) = 1/(s_i' * P_noise * s_i);
    MV_true(i) = 1/(s_i' * R_inv   * s_i);
end
for i = 1:numTest
    s_rand     = generateSparseVector(M, K);
    SM_rand(i) = 1/(s_rand' * P_noise * s_rand);
    MV_rand(i) = 1/(s_rand' * R_inv   * s_rand);
end

fprintf('MUSIC true values: [%s]\n',           num2str(SM_true',' %.2f'));
fprintf('MUSIC random:      Max = %.2f, Mean = %.2f, Median = %.2f\n', ...
        max(SM_rand), mean(SM_rand), median(SM_rand));
fprintf('MVDR true values:  [%s]\n',            num2str(MV_true',' %.2f'));
fprintf('MVDR random:       Max = %.2f, Mean = %.2f, Median = %.2f\n\n', ...
        max(MV_rand), mean(MV_rand), median(MV_rand));

% Compare peak ratios
ratioMUSIC = max(SM_true) / max(SM_rand);
ratioMVDR  = max(MV_true) / max(MV_rand);
fprintf('MUSIC true/test peak ratio = %.2f, MVDR ratio = %.2f. MUSIC better isolates the correct source vectors.\n\n', ...
        ratioMUSIC, ratioMVDR);

%% Part III: Try Again
% Reduce snapshot count to test non-invertible R
N2 = 50;
[~, ~, ~, A2]   = generateSignal(M, N2, K, PdB, Pn_dB);
[U2, Sigma2, ~] = svd(A2);
sing_vals2      = diag(Sigma2);

% Plot singular values (N=50)
figure;
stem(sing_vals2, 'filled');
title('Singular values σ_k of A  (M=100, N=50)');
xlabel('Singular-value index k');
ylabel('Singular value σ_k');
xline(3, '--', 'L=3', 'LabelOrientation','horizontal');
ylim([0, max(sing_vals2)*1.1]);
xticks(0:5:N2);
yticks(0:1:ceil(max(sing_vals2)*1.1));
grid on; box on;

fprintf('σ₃/σ₄ (N = 50) = %.3f\n\n', sing_vals2(3)/sing_vals2(4));

UL2       = U2(:,1:L);
P_noise2  = eye(M) - UL2*UL2';
SM_t2     = zeros(L,1);
SM_rand2  = zeros(L,1);

for i = 1:L
    v          = S(:,i);
    SM_t2(i)   = 1/(v' * P_noise2 * v);
end
for i = 1:numTest
    v_rand2      = generateSparseVector(M, K);
    SM_rand2(i)  = 1/(v_rand2' * P_noise2 * v_rand2);
end

fprintf('MUSIC true values (N = 50): [%s]\n',     num2str(SM_t2',' %.2f'));
fprintf('MUSIC random (N = 50):      Max = %.2f, Mean = %.2f, Median = %.2f\n\n', ...
        max(SM_rand2), mean(SM_rand2), median(SM_rand2));

fprintf('With N = 50, MUSIC still identifies the correct sources, but peaks are less pronounced and separation from random vectors is reduced.\n\n');

%% Part IV: One Last Thing
% Correlation among source vectors indicates overlap
STS = S' * S;
disp('S^T * S (source correlation matrix):');
disp(STS);

fprintf(['S^T*S diagonals = 1 (unit-norm), off-diagonals = inner products.\n', ...
         'Low off-diagonals ⇒ near-orthogonal sources (good separation);\n', ...
         'high off-diagonals ⇒ correlated sources (poorer resolution).\n']);

%% Local functions

function [S, B, V, A] = generateSignal(M, N, K, PdB, Pn_dB)
    % Build L unit-norm K-sparse source vectors
    L = numel(PdB);
    S = arrayfun(@(k) generateSparseVector(M, K), 1:L, 'UniformOutput', false);
    S = cat(2, S{:});
    % Gaussian source coefficients scaled by PdB
    B = cell2mat(arrayfun(@(p) 10^(p/20)*randn(1,N), PdB, 'UniformOutput', false)');
    % Additive white Gaussian noise
    V = 10^(Pn_dB/20) * randn(M, N);
    % Combine to form measurement matrix
    A = S*B + (1/sqrt(M))*V;
end

function vec = generateSparseVector(M, K)
    % Generate a random K-sparse unit vector
    positions = randperm(M, K);
    vec = zeros(M,1);
    vec(positions) = 1/sqrt(K);
end