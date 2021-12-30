clear all;
clc;
close all;

%% Data Generation
% no rng regulation.
fc = 2000e6;
fs = 5e9;
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
x = randi([0 3], 2000, 50);
qpskmod = comm.QPSKModulator;

for alfa = ALFA
    i = 0;
    number = number + 1;
%     for snr = SNR
        for segments = 1:2000
            i = i + 1;           
            
            waveform = qpskmod(x(segments, :)');           
            up_sig = upfirdn(waveform, raised_cos, sps);
            dn = up_sig .* carrier';
            
            em_sig = alfa{1}(1)*dn + alfa{1}(2)*dn.^2 + alfa{1}(3)*dn.^3;

            em_noisy = em_sig; %awgn(em_sig, snr, 'measured');
            
            [Bspec, waxis] = bispecd(em_noisy);
            
%             max_bs = max(max(Bspec));
%             new_bs = 1 * Bspec/max_bs;
%             flipped_bs = flip(new_bs);
            
%             folder_name = strcat('feature\', 'emitter', int2str(number), '\snr_', int2str(snr));
            folder_name = strcat('feature2\', 'emitter', int2str(number));
            mkdir(folder_name);
            
            set(gca, 'Visible', 'off');
            saveas(gcf, strcat(folder_name, '\', int2str(segments), '.png'));
%             imwrite(flipped_hs, strcat(folder_name, '\', int2str(segments), '.png'));
    
        end
%     end
end
