% Talha Akhlaq (Problem Set 5, Question 1)

clear; close all; clc;

% (a) Design each of the four filters

rp = 1.5;     % Passband ripple in dB
rs = 30;      % Stopband attenuation in dB

Fs = 40e6;             % Sampling frequency (Hz)
Fn = Fs / 2;           % Nyquist frequency (Hz)


fpass = [10e6 13e6];   % Passband: 10–13 MHz
fstop = [9e6 14e6];    % Stopband: <9 MHz or >14 MHz

% ANALOG: convert to rad/s
Wp_a = 2*pi*fpass;
Ws_a = 2*pi*fstop;

% DIGITAL: normalize to Nyquist
Wp_d = fpass / Fn;
Ws_d = fstop / Fn;

% Analog Elliptic
[n_ae, Wn_ae] = ellipord(Wp_a, Ws_a, rp, rs, 's');
[z_ae, p_ae, k_ae] = ellip(n_ae, rp, rs, Wn_ae, 'bandpass', 's');
[b_ae, a_ae] = zp2tf(z_ae, p_ae, k_ae);

% Analog Chebyshev I
[n_ac, Wn_ac] = cheb1ord(Wp_a, Ws_a, rp, rs, 's');
[z_ac, p_ac, k_ac] = cheby1(n_ac, rp, Wn_ac, 'bandpass', 's');
[b_ac, a_ac] = zp2tf(z_ac, p_ac, k_ac);

% Digital Elliptic
[n_de, Wn_de] = ellipord(Wp_d, Ws_d, rp, rs);
[z_de, p_de, k_de] = ellip(n_de, rp, rs, Wn_de, 'bandpass');
[b_de, a_de] = zp2tf(z_de, p_de, k_de);

% Digital Chebyshev I
[n_dc, Wn_dc] = cheb1ord(Wp_d, Ws_d, rp, rs);
[z_dc, p_dc, k_dc] = cheby1(n_dc, rp, Wn_dc, 'bandpass');
[b_dc, a_dc] = zp2tf(z_dc, p_dc, k_dc);

% (b) Report the orders

fprintf('\nFilter Orders:\n');
fprintf('Analog Elliptic:      Prototype = %d, Actual = %d\n', n_ae, length(p_ae));
fprintf('Analog Chebyshev I:   Prototype = %d, Actual = %d\n', n_ac, length(p_ac));
fprintf('Digital Elliptic:     Prototype = %d, Actual = %d\n', n_de, length(p_de));
fprintf('Digital Chebyshev I:  Prototype = %d, Actual = %d\n', n_dc, length(p_dc));

% (c) Plot poles and zeros for each filter

figure; zplane(z_ae, p_ae); title('Analog Elliptic – Poles & Zeros');
figure; zplane(z_ac, p_ac); title('Analog Chebyshev I – Poles & Zeros');
figure; zplane(z_de, p_de); title('Digital Elliptic – Poles & Zeros');
figure; zplane(z_dc, p_dc); title('Digital Chebyshev I – Poles & Zeros');

% (d) Plot frequency response of each filter 

% Helper for analog frequency response
function plotAnalogFR(b, a, name)
    f = linspace(0, 20e6, 1000);             % Hz
    w = 2*pi*f;                              % rad/s
    H = freqs(b, a, w);
    mag = 20*log10(abs(H));
    phase = unwrap(angle(H)) * 180/pi;
    figure;
    subplot(2,1,1); plot(f/1e6, mag); xlabel('Frequency (MHz)'); ylabel('Magnitude (dB)');
    title([name ' – Frequency Response']);
    subplot(2,1,2); plot(f/1e6, phase); xlabel('Frequency (MHz)'); ylabel('Phase (degrees)');
end

% Helper for digital frequency response
function plotDigitalFR(b, a, Fs, name)
    f = linspace(0, Fs/2, 1000);  
    [H, ~] = freqz(b, a, f, Fs);
    mag = 20*log10(abs(H));
    phase = unwrap(angle(H)) * 180/pi;
    figure;
    subplot(2,1,1); plot(f/1e6, mag); xlabel('Frequency (MHz)'); ylabel('Magnitude (dB)');
    title([name ' – Frequency Response']);
    subplot(2,1,2); plot(f/1e6, phase); xlabel('Frequency (MHz)'); ylabel('Phase (degrees)');
end

% Plot frequency responses 
plotAnalogFR(b_ae, a_ae, 'Analog Elliptic');
plotAnalogFR(b_ac, a_ac, 'Analog Chebyshev I');
plotDigitalFR(b_de, a_de, Fs, 'Digital Elliptic');
plotDigitalFR(b_dc, a_dc, Fs, 'Digital Chebyshev I');
