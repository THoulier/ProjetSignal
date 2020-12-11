clear;
close all;
clc;

%% PARTIE 2

%% Preliminaire 1

R1 = 0.98;
teta1 = 6*pi/15;

R2 = 0.98;
teta2 = 5*pi/15;

p1 = R1*exp(j*teta1);
p2 = R1*exp(-j*teta1);
p3 = R2*exp(j*teta2);
p4 = R2*exp(-j*teta2);

zeros = [1 0];
poles = [1 -p1-p2-p3-p4 p3*p4+(-p1-p2)*(-p3-p4)+p1*p2 p3*p4*(-p1-p2)+p1*p2*(-p3-p4) p1*p2*p3*p4];

N = 5000;
f = -1/2:1/N:1/2-1/N;

m = 0;
sigmacarre = 1;
bruit = sqrt(sigmacarre)*randn(1,N);
AR_process = filter(zeros, poles, bruit);

h = freqz(zeros, poles, 2*pi*f);

%% Representation temporelle

t = 0: 1/N :1-1/N;
figure;
subplot(2,1,1);
plot(t,bruit)
title('Représentation temporelle du bruit avant filtre');
xlabel('temps');
subplot(2,1,2);
plot(t,AR_process);
title('Représentation temporelle du AR');
xlabel('temps');


%% SPECTRE DE PUISSANCE
TF_AR = fftshift(fft(AR_process));
SP_AR = (1/N)*abs(TF_AR).^2;

TF_bruit = fftshift(fft(bruit));
SP_bruit = (1/N)*abs(TF_bruit).^2;


figure;
subplot(2,1,1);
plot(f,SP_AR)
title('Spectre de puissance du AR');
xlabel('frequence normalisée');
subplot(2,1,2);
plot(f,SP_bruit)
title('Spectre de puissance du bruit avant filtrage');
xlabel('frequence normalisée');

%% DENSITE SPECTRALE DE PUISSANCE
DSP_AR = SP_AR;
figure;
plot(f,DSP_AR)
title('DSP du bruit');
xlabel('frequence normalisée');
hold on;
plot(f, abs(h).^2*sigmacarre,'r','LineWidth',2)
legend('Spectre de puissance AR','DSP AR');

%% Bruit additif

eb_n0_dB = 10;
eb_n0 = 10^(eb_n0_dB/10);
Eb = sum(abs(AR_process.^2)); % Energie du filtre 
N0 = Eb/eb_n0;


m = 0;
sigmacarre = 1;
bruit_additif = sqrt(N0/2)*randn(1,N);

AR_bruite = AR_process + bruit_additif;

TF_AR_bruite = fftshift(fft(AR_bruite));
SP_AR_bruite = (1/N)*abs(TF_AR_bruite).^2;

%% Representation du SP AR et DSP AR
figure;

plot(f,SP_AR_bruite)
title('Spectre du AR bruité');
xlabel('frequence normalisée');

%% Representation temporelle de AR bruite

t = 0: 1/N :1-1/N;
figure;
subplot(2,1,1);
plot(t,AR_process)
title('Représentation temporelle du AR');
xlabel('temps');
subplot(2,1,2);
plot(t,AR_bruite);
title('Représentation temporelle du AR bruité');
xlabel('temps');
