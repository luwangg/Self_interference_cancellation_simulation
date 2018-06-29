close all;clear;clc;j=1i;
%% Parameters
N = 1e4; % Total TX symbol number
Trial = 10; % Trial number
%% Initial Parameters
H_RR_A_norm_error_average = zeros(N,1);
H_RR_AA_norm_error_average = zeros(N,1);
r_r_desired_error_average = zeros(N,1);
r_r_desired_error_rls_average = zeros(N,1);
Pre_SIR_cul = 0;
After_SIR = 0;
After_SIR_rls = 0;
%% Trial Loop
for tt = 1:Trial
%% Channel
Channel_length = 3;
H_SR = (randn(Channel_length,2)*[1;j]/sqrt(2)).'; % Rayleigh channel [1x3]
Pre_SIR = -20; % dB
H_RR = 10.^(-Pre_SIR/10/2)*(randn(Channel_length,2)*[1;j]/sqrt(2)).'; % Rayleigh channel [1x3]
%% Load Test H Data
% load('H.mat');
% Self-define
% H_SR = [0.56-j*0.43 , 0.34+j*0.33 , 0.12+j*0.22];
% H_RR = [0.45-j*0.4 , 0.32+j*0.31 , 0.1+j*0.12];
%% Noise
SNR = 20;
np = 10.^(-SNR/10); % noise power, Eb = 1
n_T = sqrt(np)*randn(N,2)*[1;j]/sqrt(2);
%% TX signal
M = 4; % QPSK
s_t = pskmod(randi([0 M-1],N,1),M,pi/4); % Un-Known signal
r_t = pskmod(randi([0 M-1],N,1),M,pi/4); % Known signal
%% RX signal
r_r = zeros(N,1);
r_r_rls = zeros(N,1);
%% Cancellation setting
A_length = 3;
A = zeros(1,A_length);          % Initial weigth
AA = zeros(1,A_length);         % Initial weigth
%% Initial Parameters
H_RR_A_norm_error = zeros(N,1);
H_RR_AA_norm_error = zeros(N,1);
desired = zeros(N,1);
interference = zeros(N,1);
%% RLS Initial
epsilon = 0.0001;
P = (epsilon^-1)*eye(A_length);
%% Main
for index = A_length:N
    desired(index) = H_SR*s_t(index:-1:index-(Channel_length-1));
    interference(index) = H_RR*r_t(index:-1:index-(Channel_length-1));
    % SGD
    r_r(index) = desired(index) + interference(index) + A*r_t(index:-1:index-(A_length-1)) + n_T(index);
    A = DF_SGD_function(A,r_r(index),r_t(index:-1:index-(A_length-1)));
    H_RR_A_norm_error(index) = 10*log10(norm(H_RR+A(1:Channel_length))^2); % dB
    % RLS
    r_r_rls(index) = desired(index) + interference(index) + AA*r_t(index:-1:index-(A_length-1)) + n_T(index);
    [AA,P] = DF_RLS_function(AA,r_r_rls(index),r_t(index:-1:index-(A_length-1)),P);
    H_RR_AA_norm_error(index) = 10*log10(norm(H_RR+AA(1:Channel_length))^2); % dB
end
H_RR_A_norm_error_average = H_RR_A_norm_error_average + H_RR_A_norm_error/Trial;
H_RR_AA_norm_error_average = H_RR_AA_norm_error_average + H_RR_AA_norm_error/Trial;
r_r_desired_error_average = r_r_desired_error_average + 10*log10(abs(r_r - desired))/Trial; % dB
r_r_desired_error_rls_average = r_r_desired_error_rls_average + 10*log10(abs(r_r_rls - desired))/Trial; % dB
Pre_SIR_cul = Pre_SIR_cul + (10*log10((norm(desired)^2) / (norm(interference)^2))) /Trial;
After_SIR = After_SIR + (10*log10((norm(desired)^2) / (norm(r_r-desired)^2))) /Trial;
After_SIR_rls = After_SIR_rls + (10*log10((norm(desired)^2) / (norm(r_r_rls-desired)^2))) /Trial;
end % Trial Loop

% Pre_SINR = 10*log10((norm(desired)^2) / ((norm(interference)^2) + norm(n_T)^2));
% After_SINR = 10*log10((norm(desired)^2) / ((norm(interference)^2) + norm(n_T)^2));
%% Plot
subplot(2,4,1);plot(H_RR_A_norm_error_average);
hold on
subplot(2,4,1);plot(H_RR_AA_norm_error_average);
hold off
title('E \{ || A[n] - H_R_R[n] ||^2 \}');xlabel('Iteration');ylabel('dB');axis square;legend('SGD','RLS');
%-----------------------------------------------------------%
subplot(2,4,2);stem(abs(H_RR));
hold on;
subplot(2,4,2);stem(abs(A(1:Channel_length)));
subplot(2,4,2);stem(abs(AA(1:Channel_length)));
hold off
title('Channel estimation');xlabel('Order');ylabel('');axis square;legend('|H_R_R|^2','|A|^2','|AA|^2');
%-----------------------------------------------------------%
subplot(2,4,3);plot(desired + interference + n_T,'.');
title('desired + interference + Noise');xlabel('');ylabel('');axis square;
%-----------------------------------------------------------%
subplot(2,4,5);plot(interference,'*');
title('interference');xlabel('');ylabel('');axis square;
%-----------------------------------------------------------%
subplot(2,4,6);plot(desired,'*');
title('desired');xlabel('');ylabel('');axis square;axis([-1.5,1.5,-1.5,1.5]);
%-----------------------------------------------------------%
subplot(2,4,7);plot(r_r,'.');
hold on
subplot(2,4,7);plot(r_r_rls,'.');
hold off
title('r_r');xlabel('');ylabel('');axis square;axis([-1.5,1.5,-1.5,1.5]);legend('SGD','RLS');
%-----------------------------------------------------------%
subplot(2,4,8);plot(r_r_desired_error_average);
hold on
subplot(2,4,8);plot(r_r_desired_error_rls_average);
hold off
title('|r_r - desired|^2');xlabel('symbol');ylabel('dB');axis square;legend('SGD','RLS');
%% Print H
% H_RR
% A
After_SIR - Pre_SIR_cul
After_SIR_rls - Pre_SIR_cul