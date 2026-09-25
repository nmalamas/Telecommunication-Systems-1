%-------- THL 1 - Ex 3 --------%
% Malamas Nikolaos(2020030180) %

clear all;close all;clc;
fgr=1;

% 1.
N = 200;
b = (sign(randn(4*N, 1)) + 1)/2;

% 2.
% Create the function bits_to_4PAM()

% 3.
Ampl = 1;
X_I = bits_to_4PAM_Ex3(b(1:2*N),Ampl)';
X_Q = bits_to_4PAM_Ex3(b(2*N+1:end),Ampl)';

% 4.
T = 0.01;
over = 10;
Ts = T/over;
Fs = 1/Ts;
hald_dur = 4;
a = 0.5;

% The time axis reffering to the two 4-PAM sequences
t_seq = 0: Ts : N*T-Ts;

% The 4-PAM sequences
XI_delta = (1/Ts) .* upsample(X_I, over);
XQ_delta = (1/Ts) .* upsample(X_Q, over);

% The SRRC filter of the two paths
[phi, t_phi] = srrc_pulse(T, over, hald_dur, a);

% By filtering the two symbol sequences 
% we essentialy execute the convolutions below
XI_t = Ts.*conv(phi,XI_delta);
XQ_t = Ts.*conv(phi,XQ_delta);

% The time axis of the convolution is given as
t_conv_start = t_phi(1) + t_seq(1);
t_conv_end = t_phi(end) + t_seq(end);
t_axis = t_conv_start : Ts : t_conv_end;

% The two waveforms in plot
figure(fgr);
subplot(2,1,1);
plot(t_axis,XI_t);
title('X_I(t)');
xlabel('time(sec)');
axis tight;
grid on;

subplot(2,1,2);
plot(t_axis,XQ_t);
title('X_Q(t)');
xlabel('time(sec)');
axis tight;
grid on;

fgr=fgr+1;

Nf = 2048;
T_total = length(XI_t)*Ts;
f_axis = [-Fs/2 : Fs/Nf : Fs/2 - 1/Nf];

% The periodogram of each waveform
FFT_XI_sqrd = abs(Ts.*fftshift(fft(XI_t,Nf))).^2;
Px_I = FFT_XI_sqrd./T_total;

FFT_XQ_sqrd = abs(Ts.*fftshift(fft(XQ_t,Nf))).^2;
Px_Q = FFT_XQ_sqrd./T_total;

figure(fgr);

subplot(2,1,1);
plot(f_axis,Px_I);
title('Periodogram of X_I');
xlabel('f (Hz)');
ylabel('P_{X,I}(f)');
grid on;
axis on;
axis tight;

subplot(2,1,2);
plot(f_axis,Px_Q);
title('Periodogram of X_Q');
xlabel('f (Hz)');
ylabel('P_{X,Q}(f)');
grid on;
axis on;
axis tight;

fgr=fgr+1;

% 5.
F0 = 200; %200 Hz

% Multiply with the corresponding carrier
X_I_mod =  2*XI_t.*cos(2*pi*F0*t_axis);
X_Q_mod = -2*XQ_t.*sin(2*pi*F0*t_axis);

% The two waveforms in plot
figure(fgr);
subplot(2,1,1);
plot(t_axis,X_I_mod);
title('X_I^{mod}(t)');
xlabel('time(sec)');
axis tight;
grid on;

subplot(2,1,2);
plot(t_axis,X_Q_mod);
title('X_Q^{mod}(t)');
xlabel('time(sec)');
axis tight;
grid on;

fgr=fgr+1;

% The periodogram of each waveform
FFT_XI_sqrd = abs(Ts.*fftshift(fft(X_I_mod,Nf))).^2;
Px_I_mod = FFT_XI_sqrd./T_total;

FFT_XQ_sqrd = abs(Ts.*fftshift(fft(X_Q_mod,Nf))).^2;
Px_Q_mod = FFT_XQ_sqrd./T_total;

figure(fgr);

subplot(2,1,1);
plot(f_axis,Px_I_mod);
title('Periodogram of X_I^{mod}');
xlabel('f (Hz)');
ylabel('P_{X,I}(f)');
grid on;
axis on;
axis tight;

subplot(2,1,2);
plot(f_axis,Px_Q_mod);
title('Periodogram of X_Q^{mod}');
xlabel('f (Hz)');
ylabel('P_{X,Q}(f)');
grid on;
axis on;
axis tight;

fgr=fgr+1;

% 6.
% The input of the channel
X_mod = X_I_mod + X_Q_mod;

% The input waveform in plot
figure(fgr);

plot(t_axis,X_mod);
title('X^{mod}(t)');
xlabel('time(sec)');
axis tight;
grid on;

fgr=fgr+1;

% The periodogram of the input waveform
FFT_XI_sqrd = abs(Ts.*fftshift(fft(X_I_mod,Nf))).^2;
Px_I_mod = FFT_XI_sqrd./T_total;

figure(fgr);
plot(f_axis,Px_I_mod);
title('Periodogram of X^{mod}');
xlabel('f (Hz)');
ylabel('P_{X,I}(f)');
grid on;
axis on;
axis tight;
fgr=fgr+1;

% 7.
% We continue by assuming that the channel
% is ideal, that is the output is the same as
% the input (impulse response δ(t))

% 8.
SNR_dB = 20;
% From the given formulas
sigma_w_sqrd = (10*Ampl^2)/(Ts*10^(SNR_dB/10));

% The WGN
mu_W = 0;
std_W = sqrt(sigma_w_sqrd);

W = std_W.*randn(1,length(X_mod)) + mu_W;

% Adding the Noise vector to the output of the channel
Y_mod = X_mod + W;

% 9.
Y_mod_I =  Y_mod.*cos(2*pi*F0*t_axis);
Y_mod_Q = -Y_mod.*sin(2*pi*F0*t_axis);

% The two waveforms in plot
figure(fgr);
subplot(2,1,1);
plot(t_axis,Y_mod_I);
title('Y_I^{mod}(t)');
xlabel('time(sec)');
axis tight;
grid on;

subplot(2,1,2);
plot(t_axis,Y_mod_Q);
title('Y_Q^{mod}(t)');
xlabel('time(sec)');
axis tight;
grid on;

fgr=fgr+1;

Nf = 2048;
T_total = length(Y_mod_I)*Ts;
f_axis = [-Fs/2 : Fs/Nf : Fs/2 - 1/Nf];

% The periodogram of each waveform
FFT_YI_sqrd = abs(Ts.*fftshift(fft(Y_mod_I,Nf))).^2;
Py_I_mod = FFT_YI_sqrd./T_total;

FFT_YQ_sqrd = abs(Ts.*fftshift(fft(Y_mod_Q,Nf))).^2;
Py_Q_mod = FFT_YQ_sqrd./T_total;

figure(fgr);

subplot(2,1,1);
plot(f_axis,Py_I_mod);
title('Periodogram of Y_I^{mod}');
xlabel('f (Hz)');
ylabel('P_{Y,I}(f)');
grid on;
axis on;
axis tight;

subplot(2,1,2);
plot(f_axis,Py_Q_mod);
title('Periodogram of Y_Q^{mod}');
xlabel('f (Hz)');
ylabel('P_{Y,Q}(f)');
grid on;
axis on;
axis tight;

fgr=fgr+1;

% 10.
% Use the same SRRC filters generated in the transmitter
[phi, t_phi_rec] = srrc_pulse(T, over, hald_dur, a);

% By filtering the two symbol sequences 
% we essentialy execute the convolutions below
YI_filtered = Ts.*conv(Y_mod_I,phi);
YQ_filtered = Ts.*conv(Y_mod_Q,phi);

% The time axis of the convolution is given as
t_conv_start = t_phi_rec(1) + t_axis(1);
t_conv_end = t_phi_rec(end) + t_axis(end);
t_axis_fltr = t_conv_start : Ts : t_conv_end;

% The two waveforms in plot
figure(fgr);
subplot(2,1,1);
plot(t_axis_fltr,YI_filtered);
title('Y_I(t)');
xlabel('time(sec)');
axis tight;
grid on;

subplot(2,1,2);
plot(t_axis_fltr,YQ_filtered);
title('Y_Q(t)');
xlabel('time(sec)');
axis tight;
grid on;

fgr=fgr+1;

Nf = 2048;
T_total = length(YI_filtered)*Ts;
f_axis = [-Fs/2 : Fs/Nf : Fs/2 - 1/Nf];

%The periodogram of each waveform
FFT_YI_sqrd = abs(Ts.*fftshift(fft(YI_filtered,Nf))).^2;
Py_I_filtered = FFT_YI_sqrd./T_total;

FFT_YQ_sqrd = abs(Ts.*fftshift(fft(YQ_filtered,Nf))).^2;
Py_Q_filtered = FFT_YQ_sqrd./T_total;

figure(fgr);
subplot(2,1,1);
plot(f_axis,Py_I_filtered);
title('Periodogram of Y_I');
xlabel('f (Hz)');
ylabel('P_{Y,I}(f)');
grid on;
axis on;
axis tight;

subplot(2,1,2);
plot(f_axis,Py_Q_filtered);
title('Periodogram of Y_Q');
xlabel('f (Hz)');
ylabel('P_{Y,Q}(f)');
grid on;
axis on;
axis tight;

fgr=fgr+1;

% 11.
% % Sample the filtered waveforms, every N-th time unit
Y_I_sampled = YI_filtered(2*hald_dur*over+1:over:end-2*hald_dur*over);
Y_Q_sampled = YQ_filtered(2*hald_dur*over+1:over:end-2*hald_dur*over);

% The filtered 16-QAM symbol sequence
symbs_recovered = [ Y_I_sampled; Y_Q_sampled ];

scatterplot(symbs_recovered');
title('The recovered symbol sequence');
xlabel('Y_{I,k}');
ylabel('Y_{Q,k}');
fgr=fgr+1;

% 12.
% Create the 'detect_4_PAM' function

% Pass into it the sample sequences to get
% the estimated sequences sent by the transmitter
est_X_I = detect_4_PAM(Y_I_sampled,hald_dur);
est_X_Q = detect_4_PAM(Y_Q_sampled,hald_dur);

% 13.
symb_err_count = 0;

for k = 1:200
    if ( (X_I(k) ~= est_X_I(k)) || (X_Q(k) ~= est_X_Q(k)) )
        symb_err_count=symb_err_count+1;
    end
end

% 14.
% Create the 'PAM_4_to_bits' function

% Unify the two sequences
est_X = [ est_X_I est_X_Q];

% Pass into it the symbol sequences to get
% the estimated bit sequence sent by the transmitter
est_bits_transmitted = PAM_4_to_bits(est_X,hald_dur);

bit_err_count = 0;

for k = 1:200
    if ( (X_I(k) ~= est_X_I(k)) || (X_Q(k) ~= est_X_Q(k)) )
        bit_err_count=bit_err_count+1;
    end
end

%%
% Symbol and Bit error probability using the Monte Carlo method
snr_max = 16;
P_Esymb = 0;
P_Ebit = 0;
P_errors = zeros(1,snr_max+2); %+2 because of snr=0

% 1.
for snr = 0 : 2 : snr_max
    symb_err_count = 0;
    bit_err_count = 0;

    for K = 1:1000
        N = 200;
        b = (sign(randn(4*N, 1)) + 1)/2;
        
        X_I = bits_to_4PAM_Ex3(b(1:2*N),hald_dur)';
        X_Q = bits_to_4PAM_Ex3(b(2*N+1:end),hald_dur)';
        
        T = 0.01;
        over = 10;
        Ts = T/over;
        Fs = 1/Ts;
        a = 0.5;
        
        XI_delta = (1/Ts) .* upsample(X_I, over);
        XQ_delta = (1/Ts) .* upsample(X_Q, over);
        
        [phi, t_phi] = srrc_pulse(T, over, hald_dur, a);
        
        XI_t = Ts.*conv(phi,XI_delta);
        XQ_t = Ts.*conv(phi,XQ_delta);
        
        F0 = 200;
        
        X_I_mod =  2*XI_t.*cos(2*pi*F0*t_axis);
        X_Q_mod = -2*XQ_t.*sin(2*pi*F0*t_axis);
        X_mod = X_I_mod + X_Q_mod;
        sigma_w_sqrd = (10*hald_dur^2)/(Ts*10^(snr/10));
        
        mu_W = 0;
        std_W = sqrt(sigma_w_sqrd);
        
        W = std_W.*randn(1,length(X_mod)) + mu_W;
        
        Y_mod = X_mod + W;
        
        Y_mod_I =  Y_mod.*cos(2*pi*F0*t_axis);
        Y_mod_Q = -Y_mod.*sin(2*pi*F0*t_axis);
        
        [phi, t_phi_rec] = srrc_pulse(T, over, hald_dur, a);
        
        YI_filtered = Ts.*conv(Y_mod_I,phi);
        YQ_filtered = Ts.*conv(Y_mod_Q,phi);
        
        Y_I_sampled = YI_filtered(2*hald_dur*over+1:over:end-2*hald_dur*over);
        Y_Q_sampled = YQ_filtered(2*hald_dur*over+1:over:end-2*hald_dur*over);
        
        symbs_recovered = [ (1/hald_dur).*Y_I_sampled; (1/hald_dur).*Y_Q_sampled ];
        
        est_X_I = detect_4_PAM(Y_I_sampled,hald_dur);
        est_X_Q = detect_4_PAM(Y_Q_sampled,hald_dur);
        for k = 1:200
            if ( (X_I(k) ~= est_X_I(k)) || (X_Q(k) ~= est_X_Q(k)) )
                symb_err_count=symb_err_count+1;
            end
        end
        
        est_X = [ est_X_I est_X_Q];
        est_bits_transmitted = PAM_4_to_bits(est_X,hald_dur);
        for k = 1:200
            if ( (X_I(k) ~= est_X_I(k)) || (X_Q(k) ~= est_X_Q(k)) )
                bit_err_count=bit_err_count+1;
            end
        end
    end
    P_Esymb = symb_err_count/(K*(length(X_I)+length(X_Q)));
    P_Ebit = bit_err_count/ (K*length(b));

    fprintf('==================================\n');
    fprintf('For SNR_dB = %d we get \n',snr);
    fprintf('Symbol Error Probability = %.3f\n',P_Esymb);
    fprintf('Bit Error Probability = %.3f\n',P_Ebit);

    % Putting the error probablities next to each other
    % symbol error at index i, correponding bit error at index i+1
    P_errors(snr+1) = P_Esymb;
    P_errors(snr+2) = P_Ebit;
end

% 2.
SNR_dB = [0:2:snr_max];
Q_arg = sqrt((1/5)*10.^(SNR_dB/10));

figure(fgr);
semilogy(SNR_dB,P_errors(1:2:end));
hold on;
semilogy(SNR_dB,3*Q(Q_arg)-(9/4)*Q(Q_arg).^2);
title('16-QAM - P_{symbol}^E');
grid on;
xlabel('SNR_{dB}');
legend('P_{symbol}^E experimental','P_{symbol}^E theoretical');
fgr=fgr+1;

% 3.

figure(fgr);
semilogy(SNR_dB,P_errors(2:2:end));
hold on;
semilogy(SNR_dB,(3/4)*Q(Q_arg));
title('16-QAM - P_{bit}^E');
grid on;
xlabel('SNR_{dB}');
legend('P_{bit}^E experimental','P_{bit}^E theoretical');
fgr=fgr+1;