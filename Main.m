close all;clear;clc;j=1i;
%% Parameters
N = 1e4; % Total TX symbol number
Trial = 10; % Trial number
%% Initial Parameters
Phi_a_SGD_norm_error_average = zeros(N,1);
Phi_a_RLS_norm_error_average = zeros(N,1);
y_desired_error_average = zeros(N,1);
y_desired_error_rls_average = zeros(N,1);
Pre_SIR_cul = 0;
After_SIR = 0;
After_SIR_rls = 0;
%% Trial Loop
for tt = 1:Trial
%% Channel
Channel_length = 3;
H = (randn(Channel_length,2)*[1;j]/sqrt(2)); % Rayleigh channel [3x1]
Pre_SIR = -20; % dB
Phi = 10.^(-Pre_SIR/10/2)*(randn(Channel_length,2)*[1;j]/sqrt(2)); % Rayleigh channel [3x1]
%% Load Test H Data
% load('H.mat');
% Self-define
% H_SR = [0.56-j*0.43 , 0.34+j*0.33 , 0.12+j*0.22];
% H_RR = [0.45-j*0.4 , 0.32+j*0.31 , 0.1+j*0.12];
%% Noise
SNR = 20;
np = 10.^(-SNR/10); % noise power, Eb = 1
n = sqrt(np)*randn(N,2)*[1;j]/sqrt(2);
%% TX signal
M = 4; % QPSK
s = pskmod(randi([0 M-1],N,1),M,pi/4); % Un-Known signal
x = pskmod(randi([0 M-1],N,1),M,pi/4); % Known signal
%% RX signal
y_SGD = zeros(N,1);
y_RLS = zeros(N,1);
%% Cancellation setting
a_length = 3;
a_GSD = zeros(a_length,1);          % Initial weigth
a_RLS = zeros(a_length,1);         % Initial weigth
%% Initial Parameters
Phi_a_SGD_norm_error = zeros(N,1);
Phi_a_RLS_norm_error = zeros(N,1);
desired = zeros(N,1);
interference = zeros(N,1);
%% RLS Initial
epsilon = 0.0001;
P = (epsilon^-1)*eye(a_length);
%% Main
for index = a_length:N
    desired(index) = H.'*s(index:-1:index-(Channel_length-1));
    interference(index) = Phi.'*x(index:-1:index-(Channel_length-1));
    % SGD
    y_SGD(index) = desired(index) + interference(index) - a_GSD.'*x(index:-1:index-(a_length-1)) + n(index);
    a_GSD = DF_SGD_function(a_GSD,y_SGD(index),x(index:-1:index-(a_length-1)));
    Phi_a_SGD_norm_error(index) = 10*log10(norm(Phi-a_GSD(1:Channel_length))^2); % dB
    % RLS
    y_RLS(index) = desired(index) + interference(index) - a_RLS.'*x(index:-1:index-(a_length-1)) + n(index);
    [a_RLS,P] = DF_RLS_function(a_RLS,y_RLS(index),x(index:-1:index-(a_length-1)),P);
    Phi_a_RLS_norm_error(index) = 10*log10(norm(Phi-a_RLS(1:Channel_length))^2); % dB
end
Phi_a_SGD_norm_error_average = Phi_a_SGD_norm_error_average + Phi_a_SGD_norm_error/Trial;
Phi_a_RLS_norm_error_average = Phi_a_RLS_norm_error_average + Phi_a_RLS_norm_error/Trial;
y_desired_error_average = y_desired_error_average + 10*log10(abs(y_SGD - desired))/Trial; % dB
y_desired_error_rls_average = y_desired_error_rls_average + 10*log10(abs(y_RLS - desired))/Trial; % dB
Pre_SIR_cul = Pre_SIR_cul + (10*log10((norm(desired)^2) / (norm(interference)^2))) /Trial;
After_SIR = After_SIR + (10*log10((norm(desired)^2) / (norm(y_SGD-desired)^2))) /Trial;
After_SIR_rls = After_SIR_rls + (10*log10((norm(desired)^2) / (norm(y_RLS-desired)^2))) /Trial;
end % Trial Loop

% Pre_SINR = 10*log10((norm(desired)^2) / ((norm(interference)^2) + norm(n_T)^2));
% After_SINR = 10*log10((norm(desired)^2) / ((norm(interference)^2) + norm(n_T)^2));
%% Plot
subplot(2,4,1);plot(Phi_a_SGD_norm_error_average);
hold on
subplot(2,4,1);plot(Phi_a_RLS_norm_error_average);
hold off
title('E \{ || a[n] - Phi[n] ||^2 \}');xlabel('Iteration');ylabel('dB');axis square;legend('SGD','RLS');
%-----------------------------------------------------------%
subplot(2,4,2);stem(abs(Phi));
hold on;
subplot(2,4,2);stem(abs(a_GSD(1:Channel_length)));
subplot(2,4,2);stem(abs(a_RLS(1:Channel_length)));
hold off
title('Channel estimation');xlabel('Order');ylabel('');axis square;legend('|Phi|^2','|a_GSD|^2','|a_RLS|^2');
%-----------------------------------------------------------%
subplot(2,4,3);plot(desired + interference + n,'.');
title('desired + interference + Noise');xlabel('');ylabel('');axis square;
%-----------------------------------------------------------%
subplot(2,4,5);plot(interference,'*');
title('interference');xlabel('');ylabel('');axis square;
%-----------------------------------------------------------%
subplot(2,4,6);plot(desired,'*');
title('desired');xlabel('');ylabel('');axis square;axis([-1.5,1.5,-1.5,1.5]);
%-----------------------------------------------------------%
subplot(2,4,7);plot(y_SGD,'.');
hold on
subplot(2,4,7);plot(y_RLS,'.');
hold off
title('r_r');xlabel('');ylabel('');axis square;axis([-1.5,1.5,-1.5,1.5]);legend('SGD','RLS');
%-----------------------------------------------------------%
subplot(2,4,8);plot(y_desired_error_average);
hold on
subplot(2,4,8);plot(y_desired_error_rls_average);
hold off
title('|r_r - desired|^2');xlabel('symbol');ylabel('dB');axis square;legend('SGD','RLS');
%% Print H
Phi
a_GSD
a_RLS
After_SIR - Pre_SIR_cul
After_SIR_rls - Pre_SIR_cul