function [b_hat] = receive(r,plot_flag)
% [b_hat] = receive(r,plot_flag)
% Receiver program for part 1 of the project. The program should produce the
% received information bits given the samples of the received signal (sampled 
% at fs Hz.)
%
% Input:
%   r = samples of the received signal
%   plot_flag = flag for plotting [0: don't plot, 1: plot] 
%
% Output:
%   b_hat = vector containing the received information bits
%
% Rev. C (VT 2016)

%********** Begin program EDIT HERE

% Complete the code below:     
m=4
%1. filter with Rx filter (Matched filter)               
disp(['Length of r before mf: ', num2str(length(r))]);
Ns = 4; % Length of the transmit pulse, defining how many samples represent each symbol
rolloffactor=0.6; % Roll-off factor for the pulse shaping filter
span=6; % Span of the filter in symbol durations
pulse = rcosdesign(rolloffactor,span,Ns(1));       % Designing a square root raised cosine filter
y = conv(r, pulse, 'same');        % Convolve the received signal with the matched filter
disp(['Length of r after mf: ', num2str(length(y))]);
%2. Sample filter output

plot_flag = 0;
fs=100000;
Ts = 1 / fs;
Rs = 1/(Ns*Ts);


y_sampled=decimate(y,4) % Puts y through a lowpass filter and keeps every n samples

disp(['Length of y_sampled after lowpass and downsample: ', num2str(length(y_sampled))]);
%3. Make decision on which symbol was transmitted
% Define your symbol amplitudes
if m==2
symbols = [-5, 5]; % Define your symbol amplitudes
 % Calculate decision boundaries (midpoints between each pair of symbols)
boundaries = (symbols(1:end-1) + symbols(2:end)) / 2;
a_hat = zeros(size(y_sampled)); % Preallocate received symbols vector
for i = 1:length(y_sampled)
    % Determine which symbol the sampled value is closest to
    [~, index] = min(abs(symbols - y_sampled(i)));
    a_hat(i) = symbols(index);    % Assign the symbol to a_hat
end
elseif m==4
    symbols=[-5,-1,1,5]
    boundaries = (symbols(1:end-1) + symbols(2:end)) / 2;
    a_hat = zeros(size(y_sampled)); % Preallocate received symbols vector
    for i = 1:length(y_sampled)
    % Determine which symbol the sampled value is closest to
    [~, index] = min(abs(symbols - y_sampled(i)));
    a_hat(i) = symbols(index);    % Assign the symbol to a_hat
    end
end
%4. Convert symbols to bits
if m==2
b_hat = zeros(size(a_hat));% Initialize the bit vector
% Define the decision boundary
decisionBoundary = 0; % Threshold between symbols to decide between bit 0 and bit 1
% Loop through received symbols and convert them to bits based on the decision boundary
for i = 1:length(a_hat)
    if a_hat(i) > decisionBoundary
        b_hat(i) = 0; % Symbol is closer to 2, so it's bit 0
    else
        b_hat(i) = 1; % Symbol is closer to 4, so it's bit 1
    end
end
elseif m == 4
     b_hat = zeros(1, 2 * length(a_hat));  % Initialize the bit vector for 4-PAM with twice the length
    % Define the decision boundaries
    % Loop through received symbols and convert them to bits based on the constellation
    for i = 1:length(a_hat)
        switch a_hat(i)
            case -5
                b_hat(2*i-1:2*i) = [0, 0];  % Map symbol -5 to bits 00
            case -1
                b_hat(2*i-1:2*i) = [0, 1];  % Map symbol -1 to bits 01
            case 1
                b_hat(2*i-1:2*i) = [1, 0];  % Map symbol 1 to bits 10
            case 5
                b_hat(2*i-1:2*i) = [1, 1];  % Map symbol 5 to bits 11
            otherwise
                error('Received symbol not in constellation');
        end
    end
end

disp(['Length of b_hat: ', num2str(length(b_hat))]); %Length of received bits after going through filters


trellis = poly2trellis(7, [171 133]); % Example trellis for a rate 1/2 code
decodedBits = vitdec(b_hat, trellis, 5, 'trunc', 'hard');






%********** DON'T EDIT FROM HERE ON
% plot Rx signals
PlotSignals(plot_flag, 'Rx', r, y, y_sampled)
%********** End program
end