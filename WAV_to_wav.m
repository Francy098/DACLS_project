% Lettura file originale (NIST SPHERE)
[data, fs] = audioread('TIMIT\TRAIN\DR1\MJEB1\SX117.WAV');

% Audiowrite genererà un file .wav
audiowrite('our_TIMIT_database\temp.wav', data, fs);