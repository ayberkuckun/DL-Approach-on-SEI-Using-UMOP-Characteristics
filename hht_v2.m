clear;
clc;
close all;

%% Data Generation
% no rng regulation.
fc = 420e3;
fs = 1e6;
t = 0:1/fs:1e-3-1/fs;
n = 1:1:457;

em1 = [1 0.5 0.3];
em2 = [1 0.08 0.6];
em3 = [1 0.01 0.01];
em4 = [1 0.01 0.4];
em5 = [1 0.6 0.08];

ALFA = {em1, em2, em3, em4, em5};

carrier = exp(1j*2*pi*fc/fs*n);

rolloff = 0.35;     % Rolloff factor
span = 8;           % Filter span in symbols
sps = 8;            % Samples per symbol

raised_cos = rcosdesign(rolloff, span, sps);

SNR = 10:2:20;
number = 0;
x = randi([0 3], 500, 50);
qpskmod = comm.QPSKModulator;

for alfa = ALFA
    i = 0;
    number = number + 1;
    for snr = SNR
        for segments = 1:500
            i = i + 1;           
            
            waveform = qpskmod(x(segments, :)');           
            up_sig = upfirdn(waveform, raised_cos, sps);
%             plot(real(up_sig)); %% pulse shape raise time?
            dn = up_sig .* carrier';

            em_sig = alfa{1}(1)*dn + alfa{1}(2)*dn.^2 + alfa{1}(3)*dn.^3;

    %         em1_sig(i) = em1(1)*dn + em1(2)*dn.^2 + em1(3)*dn.^3;
    %         em2_sig(i) = em2(1)*dn + em2(2)*dn.^2 + em2(3)*dn.^3;
    %         em3_sig(i) = em3(1)*dn + em3(2)*dn.^2 + em3(3)*dn.^3;
    %         em4_sig(i) = em4(1)*dn + em4(2)*dn.^2 + em4(3)*dn.^3;
    %         em5_sig(i) = em5(1)*dn + em5(2)*dn.^2 + em5(3)*dn.^3;

            em_noisy = awgn(em_sig, snr, 'measured');

    %         em1_noisy(i) = awgn(em1_sig(i), snr, 'measured');
    %         em2_noisy(i) = awgn(em2_sig(i), snr, 'measured');
    %         em3_noisy(i) = awgn(em3_sig(i), snr, 'measured');
    %         em4_noisy(i) = awgn(em4_sig(i), snr, 'measured');
    %         em5_noisy(i) = awgn(em5_sig(i), snr, 'measured');

            [imf, ] = emd(real(em_noisy));


    %         [imf2, ] = emd(real(em1_sig));
    %         [imf3, ] = emd(real(em1_sig));
    %         [imf4, ] = emd(real(em1_sig));
    %         [imf5, ] = emd(real(em1_sig));

            hs = hht(imf, fs, 'FrequencyResolution', 1100);
            full_hs = full(hs);
            max_hs = max(max(hs));

            new_hs = 255 * full_hs/max_hs;
            flipped_hs = uint8(floor(flip(new_hs))); %% new true way
%             map = parula;
%             imshow(flipped_hs);
%             colormap(map);
            
            folder_name = strcat('feature_x\', 'emitter', int2str(number), '\snr_', int2str(snr));
%             folder_name = strcat('feature3\', 'emitter', int2str(number));
            mkdir(folder_name);
            
            imwrite(flipped_hs, strcat(folder_name, '\', int2str(segments), '.png'));

    %         hs = hht(imf1, fs, 'FrequencyResolution', 1000);
    %         full_hs = full(hs);
    %         max_hs = max(max(hs));
    % 
    %         new_hs = 255 * full_hs/max_hs;
    %         imagesc(new_hs, [0 255]);
    %         imwrite(new_hs, "makale.png", 'Type', 'rgb');
    
        end
    end
end
