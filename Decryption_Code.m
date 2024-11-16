clc;
clear all;
close all;

% Paths to the files
enc_audio_path = "D:\24-Digital Communication Systems\Project\Audio Files\encrypted_audio-1.wav";
dec_audio_path = "D:\24-Digital Communication Systems\Project\Audio Files\decrypted_audio-1.wav";
key_path = "D:\24-Digital Communication Systems\Project\Audio Files\encryption_key-1.mat";

% Read encrypted audio and load key
[enc_audio, Fs] = audioread(enc_audio_path);
load(key_path, 'key_comb');

% Parameters for FM demodulation
carrier_freq = 1000000;  % Carrier frequency used in the modulation
kf = 200;  % Frequency deviation factor used in the modulation
t = (0:length(enc_audio)-1) / Fs;  % Time vector for demodulation

% Optional: Downsample the signal to reduce processing load
% Decimate the signal to a lower sampling rate to speed up demodulation
downsample_factor = 10;  % Adjust based on your needs (e.g., 10x downsampling)
enc_audio_downsampled = downsample(enc_audio, downsample_factor);
Fs_downsampled = Fs / downsample_factor;  % New sample rate after downsampling
t_downsampled = (0:length(enc_audio_downsampled)-1) / Fs_downsampled;

% Demodulation (FM demodulation)
demodulated_signal_downsampled = enc_audio_downsampled .* cos(2 * pi * carrier_freq * t_downsampled);

% Filter the demodulated signal to remove unwanted frequencies
[b, a] = butter(8, 0.07, 'low');  % Low-pass filter to remove high-frequency components
filt_demodulated_signal = filter(b, a, demodulated_signal_downsampled);

% Apply high-pass filter to further clean the signal (if needed)
[b, a] = butter(8, 0.02, 'high');  % High-pass filter to remove low-frequency noise
filt_demodulated_signal2 = filter(b, a, filt_demodulated_signal);

% XOR decryption with the key
dec_sig_xor = double(bitxor(int8(filt_demodulated_signal2 * 255), int8(key_comb))) / 255;

% Reverse the frequency modulation applied during encryption
dec_signal = dec_sig_xor .* cos(2 * pi * key_comb / 255);

% Normalize the decrypted signal
dec_signal = dec_signal / max(abs(dec_signal) + eps);

% If you downsampled, upsample the signal back to the original rate
if downsample_factor > 1
    dec_signal = upsample(dec_signal, downsample_factor);
end

% Save the final decrypted audio signal
audiowrite(dec_audio_path, dec_signal, Fs);

% Plot the signals
figure;
subplot(3, 1, 1);
plot(enc_audio);
title('Received Encrypted Audio');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(3, 1, 2);
plot(dec_signal);
title('Decrypted Audio');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

subplot(3, 1, 3);
plot(filt_demodulated_signal2);
title('Filtered Demodulated Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

% Play the final decrypted audio
sound(dec_signal, Fs);
disp('Playing final decrypted audio...');
pause(length(dec_signal) / Fs + 2);
