% ----THL 1 Exercise 2----- %
% Nikolaos Malamas (2020030180)

clear all;
close all;
clc
fgr=1;

%%
% A
T = 10^(-3);
over = 10;
Ts = T / over;
Fs = 1 / Ts;
A = 4;
a = 0.5;
Nf = 2048;

[phi, t_phi] = srrc_pulse(T, over, A, a);

% A.1
% Calculate the norm of FFT of the phi pulse
fft_shifted = Ts.*fftshift(fft(phi,Nf));
fft_norm = norm(fft_shifted);
fprintf('The norm of the FFT: %f \n',fft_norm);

% Draw in semilogy the ESD of the phi pulse
fft_esd = abs(fft_shifted).^2;
f_axis = [-Fs/2 : Fs/Nf : Fs/2 - 1/Nf];

figure(fgr);
semilogy(f_axis,fft_esd);
title('A1: ESD |\Phi(F)|^2');
xlabel('F(Hz)');
ylabel('|\Phi(F)|^2');
grid on;
axis on;
axis tight;
fgr=fgr+1;

% A.2
% Create the sequence of i.i.d and equiprobable bits
N = 100;
b = (sign(randn(N, 1)) + 1)/2;

% Map them using 0 --> +1 and 1 --> -1
X_n = bits_to_2PAM(b);

% Create the waveform of symbols X_n
X_delta = (1/Ts) .* upsample(X_n, over);
t_axis = 0: Ts : N*T - Ts;

X_A2 = Ts.*conv(phi,X_delta);

t_conv_start = t_phi(1) + t_axis(1);
t_conv_end = t_phi(end) + t_axis(end);
% Calculate the correct time axis for the cnovolution
t_X = [t_conv_start : Ts : t_conv_end];

figure(fgr);
plot(t_X,X_A2);
title('A2: X(t) = $\sum_{n=0}^{N-1}X_n\phi(t-nT)$', Interpreter='latex');
xlabel('t(sec)');
ylabel('X(t)');
axis tight;
axis on;
grid on;
fgr=fgr+1;

% Calculate the theoretical PSD, using the given formula
T_total = length(X_A2)*Ts;
std_dev_XA2 = mean(X_A2.^2) - mean(X_A2).^2;
Sx_theor_A2 = ( std_dev_XA2 ./T_total ).*fft_esd;

% A.3
% Calculate the FFT of sequence X_A2
FFT_X_sqrd = abs(Ts.*fftshift(fft(X_A2,Nf))).^2;

% Calculate the Periodogram of a realization of X_A2
Px = FFT_X_sqrd./T_total;

% The periodogram of a realization of X drawn using:
% Plot
figure(fgr);
plot(f_axis,Px);
title('A3: Peridogram of a realization of X in plot');
xlabel('F(Hz)');
ylabel('P_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

% Semilogy
figure(fgr);
semilogy(f_axis,Px);
title('A3: Peridogram of a realization of X in semilogy');
xlabel('F(Hz)');
ylabel('P_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;


% Redoing all the above to get 3 more realizations of X
for i=1:3
    % Create the sequence of i.i.d and equiprobable bits
    N = 100;
    b = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n = bits_to_2PAM(b);
    
    % Create the waveform of symbols X_n
    X_delta = (1/Ts) .* upsample(X_n, over);
    X = Ts.*conv(phi,X_delta);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd = abs(Ts.*fftshift(fft(X,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px = FFT_X_sqrd./T_total;
    
    % The periodogram of a realization of X drawn using:
    figure(fgr);
    
    % Plot
    subplot(2,1,1);
    plot(f_axis,Px);
    title('A3: Another realization of X using plot');
    xlabel('F(Hz)');
    ylabel('P_X(F)');
    grid on;
    axis on;
    axis tight;
    hold on;
    
    % Semilogy
    subplot(2,1,2);
    semilogy(f_axis,Px);
    title('A3: Another realization of X using semilogy');
    xlabel('F(Hz)');
    ylabel('P_X(F)');
    grid on;
    axis on;
    axis tight;
    
    fgr=fgr+1;
end

%===========================================%

% Calculate the mean of 500 Periodgrams to get PSD
Sx_exp = 0;
K = 500;
N = 100;
for k=1:K
    
    % Create the random bit sequence
    b_cur = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n_cur = bits_to_2PAM(b_cur);
    
    % Create the waveform of symbols X_n
    X_delta_cur = (1/Ts) .* upsample(X_n_cur, over);
    X_cur = Ts.*conv(phi,X_delta_cur);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd_cur = abs(Ts.*fftshift(fft(X_cur,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px_cur = FFT_X_sqrd_cur./T_total;

    Sx_exp = Sx_exp + Px_cur;
end
Sx_exp_A3 = (1/K).*Sx_exp;

figure(fgr);
semilogy(f_axis,Sx_theor_A2,'r');
hold on;
semilogy(f_axis,Sx_exp_A3,'b');
title('Experimental and Theoretical PSD in common semilogy');
legend('S_{X,theoretical}','S_{X,experimental}');
xlabel('F(Hz)');
ylabel('S_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

%===========================================%
% Redo for bigger values of K and N
% Calculate the mean of 1000 Periodgrams to get PSD
Sx_exp = 0;
K = 1000;
N = 200;
for k=1:K
    
    % Create the random bit sequence
    b_cur = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n_cur = bits_to_2PAM(b_cur);
    
    % Create the waveform of symbols X_n
    X_delta_cur = (1/Ts) .* upsample(X_n_cur, over);
    X_cur = Ts.*conv(phi,X_delta_cur);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd_cur = abs(Ts.*fftshift(fft(X_cur,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px_cur = FFT_X_sqrd_cur./T_total;

    Sx_exp = Sx_exp + Px_cur;
end
Sx_exp_A3_K1000_N200 = (1/K).*Sx_exp;

figure(fgr);
semilogy(f_axis,Sx_theor_A2,'r');
hold on;
semilogy(f_axis,Sx_exp_A3_K1000_N200,'b');
title('Experimental and Theoretical PSD in common semilogy');
legend('S_{X,theoretical}','S_{X,experimental}');
xlabel('F(Hz)');
ylabel('S_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

%===========================================%
% Redo for bigger values of K and N
% Calculate the mean of 2000 Periodgrams to get PSD
Sx_exp = 0;
K = 2000;
N = 400;
for k=1:K
    
    % Create the random bit sequence
    b_cur = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n_cur = bits_to_2PAM(b_cur);
    
    % Create the waveform of symbols X_n
    X_delta_cur = (1/Ts) .* upsample(X_n_cur, over);
    X_cur = Ts.*conv(phi,X_delta_cur);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd_cur = abs(Ts.*fftshift(fft(X_cur,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px_cur = FFT_X_sqrd_cur./T_total;

    Sx_exp = Sx_exp + Px_cur;
end
Sx_exp_A3_K2000_N400 = (1/K).*Sx_exp;

figure(fgr);
semilogy(f_axis,Sx_theor_A2,'r');
hold on;
semilogy(f_axis,Sx_exp_A3_K2000_N400,'b');
title('Experimental and Theoretical PSD in common semilogy');
legend('S_{X,theoretical}','S_{X,experimental}');
xlabel('F(Hz)');
ylabel('S_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;


% A.4
% Create the sequence of i.i.d and equiprobable bits
N = 100;
b1 = (sign(randn((N/2), 1)) + 1)/2;
b2 = (sign(randn((N/2), 1)) + 1)/2;

% Map them using 0 --> +1 and 1 --> -1
X_n = bits_to_4PAM(b1,b2);

%===========================================%

% Create the waveform of symbols X_n
X_delta = (1/Ts) .* upsample(X_n, over);
t_axis = 0: Ts : (N/2)*T - Ts;

X_A4 = Ts.*conv(X_delta,phi);

t_conv_start = t_phi(1) + t_axis(1);
t_conv_end = t_phi(end) + t_axis(end);
% Calculate the correct time axis for the cnovolution
t_X = [t_conv_start : Ts : t_conv_end];

figure(fgr);
plot(t_X,X_A4);
title('X(t) = $\sum_{n=0}^{N/2-1}X_n\phi(t-nT)$', Interpreter='latex');
xlabel('t(sec)');
ylabel('X(t)');
axis tight;
axis on;
grid on;
fgr=fgr+1;

% Calculate the theoretical PSD, using the given formula
T_total = length(X_A4)*Ts;
std_dev_XA4 = mean(X_A4.^2) - mean(X_A4).^2;
Sx_theor_A4 = ( std_dev_XA4 ./T_total ).*fft_esd;

%===========================================%

% Draw Periodogram for a realization of the new X
% Calculate the FFT of sequence X
FFT_X_sqrd = abs(Ts.*fftshift(fft(X_A4,Nf))).^2;

% Calculate the Periodogram of a realization of the new X
Px = FFT_X_sqrd./T_total;

% The periodogram of a realization of the new X drawn using:
% Plot
figure(fgr);
subplot(2,1,1);
plot(f_axis,Px);
title('Peridogram of a realization of X in plot');
xlabel('F(Hz)');
ylabel('P_X(F)');
grid on;
axis on;
axis tight;

% Semilogy
subplot(2,1,2);
semilogy(f_axis,Px);
title('Peridogram of a realization of X in semilogy');
xlabel('F(Hz)');
ylabel('P_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

% Calculate the mean of 500 Periodgrams to get PSD of the new X
Sx_exp = 0;
K = 500;
N = 100;
for k=1:K
    
    % Create the random bit sequence
    b_cur1 = (sign(randn(N, 1)) + 1)/2;
    b_cur2 = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n_cur = bits_to_4PAM(b_cur1,b_cur2);
    
    % Create the waveform of symbols X_n
    X_delta_cur = (1/Ts) .* upsample(X_n_cur, over);
    X_cur = Ts.*conv(X_delta_cur,phi);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd_cur = abs(Ts.*fftshift(fft(X_cur,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px_cur = FFT_X_sqrd_cur./T_total;

    Sx_exp = Sx_exp + (1/K).*Px_cur;
end
Sx_exp_A4 = Sx_exp;

figure(fgr);
semilogy(f_axis,Sx_theor_A4,'r');
hold on;
semilogy(f_axis,Sx_exp_A4,'b');
title('Experimental and Theoretical PSD in common semilogy');
legend('S_{X,theoretical}','S_{X,experimental}');
xlabel('F(Hz)');
ylabel('S_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

%===========================================%

figure(fgr);
semilogy(f_axis,Sx_exp_A4,'b');
hold on;
semilogy(f_axis,Sx_exp_A3,'m');
title('Comparison of the experimental PSDs');
xlabel('F(Hz)');
ylabel('S_X(F)');
legend('4-PAM S_X(F)','2-PAM: S_X(F)');
axis tight;
grid on;
fgr=fgr+1;

%===========================================%

% A.5
T_prime = 2*T;
over_prime = 2*over;
Ts = T_prime / over_prime;
Fs = 1 / Ts;
A = 4;
a = 0.5;
Nf = 2048;

% Redoing A.3, using the new T, over values
% Calculate the FFT of sequence X
FFT_X_sqrd = abs(Ts.*fftshift(fft(X_A2,Nf))).^2;

% Calculate the Periodogram of a realization of X
Px = FFT_X_sqrd./T_total;

% The periodogram of a realization of X drawn using:
% Plot
figure(fgr);
plot(f_axis,Px);
title('Peridogram of a realization of X in plot, T_{new} = 2T');
xlabel('F(Hz)');
ylabel('P_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

% Semilogy
figure(fgr);
semilogy(f_axis,Px);
title('Peridogram of a realization of X in semilogy, T_{new} = 2T');
xlabel('F(Hz)');
ylabel('P_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

% Redoing all the above to get 3 more realizations of X
for i=1:3
    % Create the sequence of i.i.d and equiprobable bits
    N = 100;
    b = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n = bits_to_2PAM(b);
    
    % Create the waveform of symbols X_n
    X_delta = (1/Ts) .* upsample(X_n, over);
    X = Ts.*conv(X_delta,phi);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd = abs(Ts.*fftshift(fft(X_A2,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px = FFT_X_sqrd./T_total;
    
    % The periodogram of a realization of X drawn using:
    figure(fgr);
    
    % Plot
    subplot(2,1,1);
    plot(f_axis,Px);
    title('Another realization of X using plot');
    xlabel('F(Hz)');
    ylabel('P_X(F)');
    grid on;
    axis on;
    axis tight;
    hold on;
    
    % Semilogy
    subplot(2,1,2);
    semilogy(f_axis,Px);
    title('Another realization of X using semilogy');
    xlabel('F(Hz)');
    ylabel('P_X(F)');
    grid on;
    axis on;
    axis tight;
    
    fgr=fgr+1;
end

%===========================================%

% Calculate the mean of 500 Periodgrams to get PSD
Sx_exp = 0;
K = 500;
N = 100;
for k=1:K
    
    % Create the random bit sequence
    b_cur = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n_cur = bits_to_2PAM(b_cur);
    
    % Create the waveform of symbols X_n
    X_delta_cur = (1/Ts) .* upsample(X_n_cur, over);
    X_cur = Ts.*conv(X_delta_cur,phi);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd_cur = abs(Ts.*fftshift(fft(X_cur,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px_cur = FFT_X_sqrd_cur./T_total;

    Sx_exp = Sx_exp + (1/K).*Px_cur;
end
Sx_exp_A5 = Sx_exp;

figure(fgr);
semilogy(f_axis,Sx_theor_A2,'r');
hold on;
semilogy(f_axis,Sx_exp_A5,'b');
title('Experimental and Theoretical PSD in common semilogy, T_{new} = 2T');
legend('S_{X,theoretical}','S_{X,experimental}');
xlabel('F(Hz)');
ylabel('S_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

%===========================================%
% Redo for bigger values of K and N
% Calculate the mean of 1000 Periodgrams to get PSD
Sx_exp = 0;
K = 1000;
N = 200;
for k=1:K
    
    % Create the random bit sequence
    b_cur = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n_cur = bits_to_2PAM(b_cur);
    
    % Create the waveform of symbols X_n
    X_delta_cur = (1/Ts) .* upsample(X_n_cur, over);
    X_cur = Ts.*conv(X_delta_cur,phi);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd_cur = abs(Ts.*fftshift(fft(X_cur,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px_cur = FFT_X_sqrd_cur./T_total;

    Sx_exp = Sx_exp + (1/K).*Px_cur;
end
Sx_exp_A5 = Sx_exp;

figure(fgr);
semilogy(f_axis,Sx_theor_A2,'r');
hold on;
semilogy(f_axis,Sx_exp_A5,'b');
title('Experimental and Theoretical PSD in common semilogy, T_{new} = 2T');
legend('S_{X,theoretical}','S_{X,experimental}');
xlabel('F(Hz)');
ylabel('S_X(F)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

%===========================================%
% Redo for bigger values of K and N
% Calculate the mean of 2000 Periodgrams to get PSD
Sx_exp = 0;
K = 2000;
N = 400;
for k=1:K
    
    % Create the random bit sequence
    b_cur = (sign(randn(N, 1)) + 1)/2;
    
    % Map them using 0 --> +1 and 1 --> -1
    X_n_cur = bits_to_2PAM(b_cur);
    
    % Create the waveform of symbols X_n
    X_delta_cur = (1/Ts) .* upsample(X_n_cur, over);
    X_cur = Ts.*conv(X_delta_cur,phi);
    
    % Calculate the FFT of sequence X
    FFT_X_sqrd_cur = abs(Ts.*fftshift(fft(X_cur,Nf))).^2;
    
    % Calculate the Periodogram of a realization of X
    Px_cur = FFT_X_sqrd_cur./T_total;

    Sx_exp = Sx_exp + (1/K).*Px_cur;
end
Sx_exp_A5 = Sx_exp;

figure(fgr);
semilogy(f_axis,Sx_theor_A2,'r');
hold on;
semilogy(f_axis,Sx_exp_A5,'b');
title('Experimental and Theoretical PSD in common semilogy, T_{new} = 2T');
legend('S_{X,theoretical}','S_{X,experimental}');
ylabel('S_X(F)');
xlabel('F(Hz)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

% B
% Create 5 realizations of Y,
% by combining different values of X and Phi

% Define t axis
step = 0.01;
t = [-3*2*pi : step : 3*2*pi];

% Define the frequency
F0 = 1/(2*pi);

% Y1
X = 5;
Phi = 0;
Y1 = X*cos(2*pi*F0*t + Phi);

% Y2
X = 5;
Phi = pi/2;
Y2 = X*cos(2*pi*F0*t + Phi);

% Y3
X = 10;
Phi = 0;
Y3 = X*cos(2*pi*F0*t + Phi);

% Y4
X = 10;
Phi = pi/2;
Y4 = X*cos(2*pi*F0*t + Phi);

% Y5
X = 20;
Phi = pi;
Y5 = X*cos(2*pi*F0*t + Phi);

figure(fgr);

plot(t,Y1);
hold on;
plot(t,Y2);
hold on;
plot(t,Y3);
hold on;
plot(t,Y4);
hold on;
plot(t,Y5);

title('5 Realizations of Y');
legend('Y(t) = 5cos(2\piF_0t)','Y(t) = 5cos(2\piF_0t+\pi/2)','Y(t) = 10cos(2\piF_0t)','Y(t) = 10cos(2\piF_0t+\pi/2)','Y(t) = 20cos(2\piF_0t+\pi)');
xlabel('t(sec)');
ylabel('Y(t)');
axis tight;
grid on;
axis on;
fgr=fgr+1;
