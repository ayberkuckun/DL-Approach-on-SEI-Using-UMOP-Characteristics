clear all;
clc;
close all;

%% Data Generation
rng('default');

% Fs = 1e8;
% snr = -10:1:20;
% pw = 2.6e-6;
% t = 0:1/Fs:3*pw;

% xt = (A + dA).*exp(1j*(2*pi*fc*t));
% yt = awgn(xt, snr(31));

%% Creating sinusoidal Pulse and its FFT Spectogram
 
A = 1;
fc = 710e6;
freq_sim = 1e10;    % Simulation Frequency
ts = 1/freq_sim;   % Period
PW = 2.6e-6;         % Pulsewidth (1)
%PW = 3e-6;        % Pulsewidth (2)
%PW = 5e-7;        % Pulsewidth (3)
t = ts:ts:PW;      % Time
mt = cos(2*pi*5000*t);
kp = 100;
tetha = kp*mt;
% N = 1e6;
% freq_pulse = 710e6; % Pulse Frequency

dA = (rand(size(t))*2-1) * 0.005; % linspace(-0.9, 0.9, 10) % 0.9*cos(2*pi*fa*t + ka*pi*t^2)

dtetha = (rand(size(t))*2-1) * pi/100; % linspace(-pi/36, pi/36, 10);

% xt = (A).*exp(1j*(2*pi*fc*t));
% xt = (A + dA).*exp(1j* dtetha);
xt = (A + dA).*exp(1j*(2*pi*fc*t + tetha + dtetha));

trap = 0:1/3000:1-1/3000;
trap2 = ones(1, 10000*2);
trap3 = flip(trap);
trapezoid = [trap trap2 trap3];

dump = 0.8;
models = 1*(1- 1/sqrt(1-dump^2)*exp(-dump*2*pi*fc*t).* sin(2*pi*fc*sqrt(1-dump^2)*t + atan(sqrt(1-dump^2)/dump)));
dene2 = trapezoid .* models;
plot(dene2);


dene = xt .* trapezoid;
plot(t, dene);
% plot(trapezoid);

% damped = 

% trap_f = 1*PW*abs(sin(pi*PW*))*abs()

% zt = (A + dA).*exp(1j*(2*pi*fc*t + dtetha));
plot(real(xt));
yt = awgn(dene, 25, 'measured');
% zyt = awgn(zt, 25, 'measured');
noise = awgn(zeros(size(real(xt))), 25);

t_new = ts:ts:PW*3;

figure(1);
pulse = [noise yt noise];  % Pulse
plot(t_new, real(pulse));

% figure(2);
% zpulse = [noise zyt noise];  % Pulse
% plot(real(zpulse));
% % hm = hilbert(pulse);
% iam2 = abs(hm);

% plot(t_new, iam2)
% title(['Sinusoidal pulse signal PW = ' num2str(PW*1e6) 'us'])
% xlabel('Time')
% ylabel('Amplitude')


% zt = [zeros(1, 10) xt zeros(1, 10)];
% 
% figure(1);
% plot(real(zt));

figure(3);
% [upper, lower] = envelope(real(pulse));
envelope(real(pulse));
% plot(t_new, upper);

% pulsewidth(real(xt));

%% Instantaneous Amplitude and Frequency
% hx = hilbert(xt);
% iam = abs(hx);
% ifq = Fs/(2*pi)*diff(unwrap(angle(hx)));
% 
% figure(6);
% plot(t, iam);
% figure(7);
% plot(t(1:(length(t)-1)), ifq);

%% Envelope
% zt = (A + dA/10).*exp(1j*dtetha);
% figure(2);
% pw = pulsewidth(real(zt), Fs); %önce compase et
% pulse = zt(1:pw*10000);
% figure(3);
% envelope(abs?(real(zt))); %rectify??
% plot(t, upper);
% hm = hilbert(mt);
% iam2 = abs(hm);
% figure(66);
% plot(t, iam2);
% figure(663);
% pw2 = pulsewidth(mt, Fs); %önce compase et
% t0 = pw2*10000*1.5;
% t1 = pw2*10000*2.5;
% pulse2 = mt(100:300);
% envelope(pulse2);
% figure(6635);
% [upper2, ] = envelope(pulse2);
% plot(t(1:201), upper2);
% stft(pulse, Fs);