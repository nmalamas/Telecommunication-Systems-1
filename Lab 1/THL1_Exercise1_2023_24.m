% ----THL 1 Exercise 1----- %
% Nikolaos Malamas (2020030180)

clear all;
close all;
fgr=1;

% Θ.1
T = 10;
step=0.01;
tau = -2*T:step:2*T;
R=((tau)./T +1).*(-T < tau & tau <= 0) + ((-tau)./T +1).*(0 < tau & tau <= T);

figure(fgr);
plot(tau,R);
title('\bf{R_{\phi\phi}(\tau)} for \phi(t)');
legend('R_{\phi\phi}(\tau)');
xlabel ('\tau','FontSize',14);
axis on;
grid on;
fgr=fgr+1;

% Θ.2
R=((tau)./T +1).*(-T < tau & tau <= 0) + ((-tau)./T +1).*(0 < tau & tau <= T);
figure(fgr);
plot(tau,R);
title('\bf{R_{\phi\phi}(\tau)} for \phi(t-2)');
legend('R_{\phi\phi}(\tau)');
xlabel ('\tau','FontSize',14);
axis on;
grid on;
fgr=fgr+1;

% Θ.3
R = ((-tau)./T -1).*(-T < tau & tau <= -T/2) + ((3*tau./T)+1).*(-T/2 < tau & tau <= 0) ...
  + ((-3*tau./T)+1).*(0 < tau & tau <= T/2)  + ((tau)./T -1).*(T/2 < tau & tau <= T);
figure(fgr);
plot(tau,R);
title('\bf{R_{\phi\phi}(\tau)} for \phi(t)');
legend('R_{\phi\phi}(\tau)');
xlabel ('\tau','FontSize',14);
axis on;
grid on;
fgr=fgr+1;

% A.1
T=10^(-2);
over=10;
Ts=T/over;
A=4;

% i)

% Generate the pulse for a=0
a=0;
T=10^(-2);

[phi_a0, t] = srrc_pulse(T, over, A, a);

% Generate the pulse for a=0.5
a=0.5;
T=10^(-2);

[phi_a05, t] = srrc_pulse(T, over, A, a);

% Generate the pulse for a=1
a=1;
T=10^(-2);

[phi_a1, t] = srrc_pulse(T, over, A, a);

% Common plot for the phi's
figure(fgr);
plot(t, phi_a0);
hold on;
plot(t, phi_a05);
hold on;
plot(t, phi_a1);
hold on;

legend('\phi(t) for \alpha=0', '\phi(t) for \alpha=0.5', '\phi(t) for \alpha=1');
axis tight;
xlabel('t(sec)','FontSize', 12);
grid on;
axis on;
fgr=fgr+1;

% A.2
Fs=1/Ts;
Nf=1024;
F_axis=-Fs/2:Fs/Nf:Fs/2-1/Nf;

% FFT of phi for a=0
fft_phi_a0 = fftshift(Ts.*fft(phi_a0, Nf));%evala Ts*
esd_phi_a0 = abs(fft_phi_a0).^2;

% FFT of phi for a=0.5
fft_phi_a05 = fftshift(Ts.*fft(phi_a05, Nf));%evala Ts*
esd_phi_a05 = abs(fft_phi_a05).^2;

% FFT of phi for a=1
fft_phi_a1 = fftshift(Ts.*fft(phi_a1, Nf));%evala Ts*
esd_phi_a1 = abs(fft_phi_a1).^2;

% (a)
figure(fgr);
plot(F_axis, esd_phi_a0);
hold on;
plot(F_axis, esd_phi_a05);
hold on;
plot(F_axis, esd_phi_a1);
hold on;

legend('|\Phi(F)|^2 for \alpha=0', '|\Phi(F)|^2 for \alpha=0.5', '|\Phi(F)|^2 for \alpha=1');
axis tight;
xlabel('F(Hz)');
grid on;
axis on;
fgr=fgr+1;

% (b)
figure(fgr);
semilogy(F_axis, 1/norm(esd_phi_a0) .* esd_phi_a0);
hold on;
semilogy(F_axis, 1/norm(esd_phi_a05) .* esd_phi_a05);
hold on;
semilogy(F_axis, 1/norm(esd_phi_a1) .* esd_phi_a1);
hold on;

legend('|\Phi(F)|^2 for \alpha=0', '|\Phi(F)|^2 for \alpha=0.5', '|\Phi(F)|^2 for \alpha=1');
axis tight;
xlabel('F(Hz)');
grid on;
axis on;
fgr=fgr+1;

% A.3
c = T/10^3 .* ones(length(F_axis));

figure(fgr);
semilogy(F_axis, 1/norm(esd_phi_a0) .* esd_phi_a0);
hold on;
semilogy(F_axis, 1/norm(esd_phi_a05) .* esd_phi_a05);
hold on;
semilogy(F_axis, 1/norm(esd_phi_a1) .* esd_phi_a1);
hold on;
semilogy(F_axis, c, 'LineWidth', 2);
hold on;

legend('|\Phi(F)|^2 for \alpha=0', '|\Phi(F)|^2 for \alpha=0.5', '|\Phi(F)|^2 for \alpha=1');
axis tight;
axis([-500 500 10^-7 10^3]);
xlabel('F(Hz)');
grid on;
axis on;
fgr=fgr+1;

%-----------------------------------------%

c = T/10^5 .* ones(length(F_axis));

figure(fgr);
semilogy(F_axis, 1/norm(esd_phi_a0) .* esd_phi_a0);
hold on;
semilogy(F_axis, 1/norm(esd_phi_a05) .* esd_phi_a05);
hold on;
semilogy(F_axis, 1/norm(esd_phi_a1) .* esd_phi_a1);
hold on;
semilogy(F_axis, c, 'LineWidth', 2);
hold on;

legend('|\Phi(F)|^2 for \alpha=0', '|\Phi(F)|^2 for \alpha=0.5', '|\Phi(F)|^2 for \alpha=1');
axis tight;
axis([-500 500 10^-9 10^3]);
xlabel('F(Hz)');
grid on;
axis on;
fgr=fgr+1;

% B.1
%%
T = 10^(-2);
A = 4;

for a=0:0.5:1
    [phi, t] = srrc_pulse(T, over, A, a);

    for k=0:1:2*A

        figure(fgr);
        
        subplot(2,1,1);
        plot(t+k*T, phi);
        hold on;
        plot(t, phi);
        title(['\alpha= ',num2str(a), ' and k=',num2str(k)]);
        legend('\phi(t)',' \phi(t-kT)');
        xlabel('t(sec)','FontSize', 12);
        grid on;
        axis on;
        axis tight;
        
        subplot(2,1,2);
        phi_cut = [zeros(1, k*over) phi(1: end- k*over)];
        product = phi .* phi_cut;
        plot(t, product);
        legend('\phi(t)\times\phi(t-kT)');
        xlabel('t(sec)','FontSize', 12);
        grid on;
        axis on;
        axis tight;
        
        fgr=fgr+1;

        integral_approx = sum(product.*Ts);
        fprintf('For k=%d, a=%.1f the inner product of phi and phi shifted by kT is ~%f \n',k,a,integral_approx);
        if(k == 2*A)
            fprintf('-----------------------------------------\n');
        end
    end
end
%%

% C
T = 10^(-2);
over = 10;
a = 0.5;
A = 4;

% C.1
N = 100;
b = (sign(randn(N, 1)) + 1)/2;

% C.2
% (a)
X_k = bits_to_2PAM(b);

% (b)
% Simulate the 2-PAM symbol sequence
X_delta = Ts * upsample(X_k, over);

% Plot the 2-PAM symbol sequence
step = 1/over;
t_axis = 0: Ts : N*T - Ts;

figure(fgr);
plot(t_axis,X_delta);
xlabel('t (sec)','FontSize', 12);
legend('X_{\delta}(t)');
title('2-PAM symbol sequence');
fgr=fgr+1;

% (c)
T = 10^(-2);
over = 10;
a = 0.5;
A = 4;

% Create the SRRC pulse and its time axis
[phi, t_phi] = srrc_pulse(T, over, A, a);

% Simulate the convolution between phi and X_delta
X = conv(X_delta,phi);

t_conv_start = t_phi(1) + t_axis(1);
t_conv_end = t_phi(end) + t_axis(end);
% Calculate the correct time axis for the cnovolution
t_X = [t_conv_start : Ts : t_conv_end];

figure(fgr);
plot(t_X,X);
axis tight;
grid on;
xlabel('t (sec)','FontSize', 12);
title('X(t)=$ X_{\delta} * \phi(t)$','Interpreter','latex');
fgr=fgr+1;

% (d)
phi_rev = fliplr(phi);

% Simulate the convolution between X and phi_reverersed
Z = conv(X,phi_rev);

% Since phi is an SRRC, its reversed will its non-zero values
% in the same time interval. 
t_phi_rev = t_phi;

t_conv_start = t_phi_rev(1) + t_X(1);
t_conv_end = t_phi_rev(end) + t_X(end);
t_Z = [t_conv_start : Ts : t_conv_end];

%%
figure(fgr);
plot(t_Z,Z);
axis tight;
grid on;
xlabel('t (sec)','FontSize', 12);
title('Z(t)=$ X(t) * \phi(-t)$','Interpreter','latex');
fgr=fgr+1;

for k=0:N-1
    fprintf('For k=%d, Z(kT)=%.2f and Xk=%d \n',k,Z((k+1)*over),X_k(k+1));
end

%%
figure(fgr);
plot(t_Z,Z);hold on;
title('Z(t) and the symbol sequence X');
stem([0 : N-1]*T, X_k, 'filled', 'LineStyle',':',...
     'marker','o','markersize',4,'MarkerFaceColor','red');
grid on;
xlabel('t (sec)','FontSize', 12);
axis tight;