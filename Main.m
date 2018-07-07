close all;clear;clc;j=1i;
%% Parameters
N = 1e3; % Total TX symbol number (Iteration number)
Trial = 1000; % Trial number
%% Initial Parameters
% Channel Norm error
Psi_a_norm_error_average_SGD = zeros(N,1);
Psi_a_norm_error_average_RLS = zeros(N,1);
% error
Desired_y_error_average_SGD = zeros(N,1);
Desired_y_error_average_RLS = zeros(N,1);
% SIR
SIR_average_SGD = zeros(N,1);
SIR_average_RLS = zeros(N,1);
%% Trial Loop
for tt = 1:Trial
    %% Signal
    M = 4; % QPSK
    s = pskmod(randi([0 M-1],N,1),M,pi/4); % Un-Known signal
    x = pskmod(randi([0 M-1],N,1),M,pi/4); % Known signal
    %% Channel
    Channel_length = 3;
    H = (randn(Channel_length,2)*[1;j]/sqrt(2)); % Rayleigh channel [3x1]
    SIR_pre = -30; % dB
    Psi = 10.^(-SIR_pre/10/2)*(randn(Channel_length,2)*[1;j]/sqrt(2)); % Rayleigh channel [3x1]
    %% Load Test H Data
    % load('H.mat');
    % Self-define
    % H_SR = [0.56-j*0.43 , 0.34+j*0.33 , 0.12+j*0.22];
    % H_RR = [0.45-j*0.4 , 0.32+j*0.31 , 0.1+j*0.12];
    %% AWGN Noise
    SNR = 20;
    np = 10.^(-SNR/20); % noise power, Eb = 1 dB
    n = sqrt(np)*randn(N,2)*[1;j]/sqrt(2);
    %% Trial Loop Initial Parameters
    Desired = zeros(N,1);
    Self_Interference = zeros(N,1);

    % Channel Norm error
    Phi_a_norm_error_SGD = zeros(N,1);
    Phi_a_norm_error_RLS = zeros(N,1);

    % SIR
    SIR_SGD = zeros(N,1);
    SIR_RLS = zeros(N,1);
    
    % Cancellation setting
    a_length = 3;
    a_GSD = zeros(a_length,1);  % Initial weigth
    a_RLS = zeros(a_length,1);  % Initial weigth
    y_SGD = zeros(N,1);
    y_RLS = zeros(N,1);

    % RLS Initial
    epsilon = 0.0001;
    P = (epsilon^-1)*eye(a_length);
    %% Main
    for index = a_length:N
        Desired(index) = H'*s(index:-1:index-(Channel_length-1));
        Self_Interference(index) = Psi'*x(index:-1:index-(Channel_length-1));
        
        % SGD
        y_SGD(index) = Desired(index) + Self_Interference(index) - a_GSD'*x(index:-1:index-(a_length-1)) + n(index);
        a_GSD = DF_SGD_function(a_GSD,y_SGD(index),x(index:-1:index-(a_length-1)));
        Phi_a_norm_error_SGD(index) = 10*log10(norm(Psi-a_GSD(1:Channel_length))^2); % dB
        
        % RLS
        y_RLS(index) = Desired(index) + Self_Interference(index) - a_RLS'*x(index:-1:index-(a_length-1)) + n(index);
        [a_RLS,P] = DF_RLS_function(a_RLS,y_RLS(index),x(index:-1:index-(a_length-1)),P);
        Phi_a_norm_error_RLS(index) = 10*log10(norm(Psi-a_RLS(1:Channel_length))^2); % dB
    end

    Psi_a_norm_error_average_SGD = Psi_a_norm_error_average_SGD + Phi_a_norm_error_SGD/Trial; % dB
    Psi_a_norm_error_average_RLS = Psi_a_norm_error_average_RLS + Phi_a_norm_error_RLS/Trial; % dB

    Desired_y_error_average_SGD = Desired_y_error_average_SGD + 20*log10(abs(y_SGD - Desired))/Trial; % dB
    Desired_y_error_average_RLS = Desired_y_error_average_RLS + 20*log10(abs(y_RLS - Desired))/Trial; % dB

    SIR_average_SGD = SIR_average_SGD + 20*log10(abs(Desired)./abs(Desired-y_SGD)) /Trial; % dB
    SIR_average_RLS = SIR_average_RLS + 20*log10(abs(Desired)./abs(Desired-y_RLS)) /Trial; % dB
end % Trial Loop
%% Plot
figure(1),plot(Psi_a_norm_error_average_SGD);
hold on
figure(1),plot(Psi_a_norm_error_average_RLS);
hold off
title('E \{ || \Psi - \alpha_n ||^2 \}');xlabel('Iteration (Symbol)');ylabel('dB');legend('SGD','RLS');
%-----------------------------------------------------------%
figure(2),plot(Desired_y_error_average_SGD);
hold on
figure(2),plot(Desired_y_error_average_RLS);
hold off
title('|d_n - y_n|^2');xlabel('Iteration (Symbol)');ylabel('dB');legend('SGD, \mu = 0.01','RLS');
%-----------------------------------------------------------%
figure(3),plot([1,N],[SNR,SNR],':','Color',[0.4940, 0.1840, 0.5560],'Linewidth',2);
hold on
figure(3),plot(SIR_average_RLS);
figure(3),plot(SIR_average_SGD,'Color',[0, 0.4470, 0.7410]);
hold off
title('SINR');xlabel('Iteration (Symbol)');ylabel('dB');legend('Self-Interference free floor','RLS','SGD, \mu = 0.01');