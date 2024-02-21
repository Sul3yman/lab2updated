function s = transmit(b,plot_flag)
% s = transmit(b,plot_flag)
% Transmitter program for part 1 of the project. The program should produce samples
% of the transmitted signal. The sample rate is fs Hz.
%
% Input:
%   b = vector containing the information bits to be transmitted
%   plot_flag = flag for plotting [0: don't plot, 1: plot]  
%
% Output:
%   s = vector containing samples of the transmitted signal at at rate of fs Hz
%
% Rev. C (VT 2016)

%********** Begin program, EDIT HERE

% Complete the code below to create samples of the transmitted signal.
m=4
%1. Convert bits to symbols
% Define the symbol rate Rs (symbols per second)
fs=100000;
Ns = 4;
Ts=1/fs
Rs = 1/(Ns*Ts); % Example: 1000 symbols per second, defining the rate at which symbols are transmitted



trellis = poly2trellis(7, [171 133]); % Example trellis for a rate 1/2 code
encodedBits = convenc(b, trellis);

if m==2
constellation = [-5 5];         % Specify constellation here (vector)
a = zeros(size(encodedBits)); % Initialize the symbol vector 'a' with the same size as bit vector 'b'
for i = 1:length(encodedBits)
    a(i) = constellation(encodedBits(i) + 1); % Convert the bits in vector b to symbols in vector a based on the constellation mapping
end
elseif m == 4
    constellation = [-5, -1, 1, 5];  % 4-PAM constellation
    a = zeros(1, length(encodedBits)/2);  % Since each symbol represents two bits, the size is half
    % Convert pairs of bits in vector b to symbols in vector a
    for i = 1:2:length(encodedBits)-1
        % Map the pair of bits to an integer index: 00->1, 01->2, 10->3, 11->4
        index = 2*encodedBits(i) + encodedBits(i+1) + 1;  
        a((i+1)/2) = constellation(index);
    end
end
      % Convert the bits in vector b to symbols in vector a
plot_flag = 0;
%2. Pulse Amplitude Modulation

rolloffactor=0.6; % Define the roll-off factor for the pulse shaping filter
span=6; % Define the span of the pulse shaping filter in symbol durations
pulse = rcosdesign(rolloffactor,span,Ns(1),'sqrt');% Specify the transmit pulse here (vector)



s = zeros(1,length(a)*Ns);    % Initialize the transmit signal vector 's' to accommodate all samples
for i = 1:length(a)          %Calculate the start index for the current symbol in the transmit signal vector
    startIndex=(i-1)*Ns(1)+1 %Calculate the end index for the current symbol in the transmit signal vector
    endIndex=i*Ns(1)
    s(startIndex:endIndex) = a(i) * pulse(1:Ns(1)); %  Modulate the current symbol onto the pulse and insert it into the transmit signal vector
end
disp(['Length of s: ', num2str(length(s))]); % Length of the transmitted signal
 

disp(['Length of b: ', num2str(length(encodedBits))]); %Length of the transmitted bits
%********** DON'T EDIT FROM HERE ON
% plot Tx signals
PlotSignals(plot_flag, 'Tx', a, s)
%********** End program
end