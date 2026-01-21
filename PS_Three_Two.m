%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROBLEM 2
% 
% We represent a finite-duration discrete-time vector with:
%   - x: a row vector of signal values
%   - n0: the time index for the first entry in x
%
% Then we write two functions:
%   (a) plotSignal(sig)
%       - Creates a stem plot with the time axis labeled properly
%       - Extends the time axis two points to the left and right
%       - Plots 0 outside the signal's support
%
%   (b) convSignal(sig1, sig2)
%       - Uses MATLAB conv function to compute the convolution
%       - Determines the correct start time for the result
%
%   (c) We apply these functions using h, x from the previous problem
%       and create stem plots for h, x, y = h * x as separate subplots.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (c) MAIN SCRIPT: Define h, x, compute y, and plot all three.

% From the previous problem, the signals are:
%   h = [2 1 4 5 2] with underscore at n=0
%   x = [4 2 3 2 1] with underscore at n=0

% Define h
h.x  = [2 1 4 5 2];
h.n0 = 0;  % first sample corresponds to n=0

% Define x
x.x  = [4 2 3 2 1];
x.n0 = 0;  % first sample corresponds to n=0

% Compute convolution: y = h * x
y = convSignal(h, x);

% Create a figure with subplots for h, x, and y
figure;

subplot(3,1,1);
plotSignal(h);
title('Signal h');

subplot(3,1,2);
plotSignal(x);
title('Signal x');

subplot(3,1,3);
plotSignal(y);
title('Convolution y = h * x');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (a) FUNCTION: plotSignal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotSignal(sig)
    % plotSignal(sig)
    % Creates a stem plot for a finite-duration discrete-time signal.
    %
    %   sig.x  -> row vector of signal values
    %   sig.n0 -> integer index for the first sample in sig.x
    %
    % The time axis is extended two points to the left and right,
    % and zeros are shown outside the signal's support.

    x  = sig.x;
    n0 = sig.n0;
    N  = length(x);

    % Define extended time axis: from (n0 - 2) to (n0 + N - 1 + 2)
    n = (n0 - 2):(n0 + N - 1 + 2);

    % Initialize all values to zero
    y = zeros(1, length(n));

    % The index corresponding to n0 is always 3 in this extended vector:
    %   Because n0 - (n0 - 2) = 2, then +1 => index = 3.
    startIndex = 3;

    % Place the original signal samples in the correct positions
    y(startIndex : startIndex + N - 1) = x;

    % Stem plot
    stem(n, y, 'filled');
    xlabel('Time index n');
    ylabel('Amplitude');
    grid on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (b) FUNCTION: convSignal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = convSignal(sig1, sig2)
    % z = convSignal(sig1, sig2)
    % Convolves two finite-duration signals (sig1, sig2).
    %
    %   sig1.x, sig1.n0
    %   sig2.x, sig2.n0
    %
    % The result z has:
    %   z.x  = conv(sig1.x, sig2.x)
    %   z.n0 = sig1.n0 + sig2.n0  (start index for the convolution)

    z.x = conv(sig1.x, sig2.x);
    z.n0 = sig1.n0 + sig2.n0;
end
