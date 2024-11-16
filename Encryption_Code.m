clc;
clear all;
close all;

% Paths to the audio files
in_audio = "D:\24-Digital Communication Systems\Project\Audio Files\sample_audio-1.wav";
enc_audio = "D:\24-Digital Communication Systems\Project\Audio Files\encrypted_audio-1.wav";
key_file = "D:\24-Digital Communication Systems\Project\Audio Files\encryption_key-1.mat";

% Read the input audio
[audio, Fs] = audioread(in_audio);
len = length(audio);

% Parameters for chaotic sequences
r1 = 3.99; 
r2 = 3.8;
r3 = 3.1;

% Initialize key components
k1 = zeros(len, 1);
k2 = zeros(len, 1);
k3 = zeros(len, 1);

% Set initial values for chaotic sequences
k1(1) = 0.5;
k2(1) = 0.3;
k3(1) = 0.9;

% Generate chaotic key components for each sequence
for i = 2:len
    k1(i) = r1 * k1(i-1) * (1 - k1(i-1));
    k2(i) = r2 * k2(i-1) * (1 - k2(i-1));
    k3(i) = r3 * k3(i-1) * (1 - k3(i-1));
end

% Normalize and scale the chaotic keys
k1 = round((k1 - min(k1)) * 255);
k2 = round((k2 - min(k2)) * 255);
k3 = round((k3 - min(k3)) * 255);

% Combine the chaotic key components into one key
key_comb = k1 + k2 + k3;

% FM Encryption: Encrypt the audio using the chaotic keys
fm_encrypted = audio .* cos(2 * pi * key_comb / 255);  % Apply FM modulation to the audio

% XOR encryption with the key_comb
enc_sig = bitxor(int8(fm_encrypted * 255), int8(key_comb));
enc_sig = double(enc_sig) / 255;

% FM Modulation after encryption (Apply FM to the encrypted signal)
% Carrier frequency for FM modulation (choose a high enough carrier frequency)
carrier_freq = 1000000;  % 1 MHz carrier frequency for FM
kf = 2 * pi * 50000;  % Frequency deviation constant (adjust as needed)

% Create time vector based on the length of the encrypted signal
t = (0:length(enc_sig)-1)' / Fs;  % Ensure t is a column vector

% Apply FM modulation to the encrypted signal
fm_modulated_signal = cos(2 * pi * carrier_freq * t + kf * enc_sig);  % Ensure matching dimensions

% Add AWGN noise to the FM modulated signal
SNR = 20;  % Signal-to-noise ratio
trans_sig = awgn(fm_modulated_signal, SNR, 'measured');

% Normalize the transmitted signal to be within the valid range for audio
trans_sig = trans_sig / max(abs(trans_sig));  % Normalize the signal to [-1, 1]

% Save the encrypted audio and key
audiowrite(enc_audio, trans_sig, Fs);  % Now the data is within the valid range

% Play the original audio
sound(audio, Fs);
disp('(Original Audio) Playing...');
pause(length(audio) / Fs + 2);

% Play the transmitted audio (with noise)
sound(trans_sig, Fs);
disp('(Transmitted Audio with Noise) Playing...');
pause(length(trans_sig) / Fs + 2);

% Plot the signals
figure;
subplot(4, 1, 1);
plot(audio);
title('Original Audio');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(4, 1, 2);
plot(fm_encrypted);
title('FM Encrypted Audio');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(4, 1, 3);
plot(fm_modulated_signal);
title('FM Modulated Encrypted Audio');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(4, 1, 4);
plot(trans_sig);
title('Transmitted Audio (with AWGN)');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
