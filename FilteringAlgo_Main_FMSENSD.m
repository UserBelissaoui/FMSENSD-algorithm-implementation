%% Implementation of Minimal Mean Squared Error in Frequency domain with N-coefficient and Signal Detection (FMSENSD) Algorithm
%** MATLAB version: R2019b
%** Abderrahim BELISSAOUI
% Date :    2/4/2021                 
%% Clear command window, close plots ...
clear all;  %#ok<*CLALL>
close all;
clc ;
%% Parameters and audio file reading
%% Sounds signals

    %* x : the mixed sound signal
    %* r : the reference noise signal (background signal)
    %* e : output signal  (reconstructed signal)
[x, Fs_x] = audioread('D:\Noise filtering using Adaptive algorithm\Dataset\mixed1.wav'); %  
[r, Fs_r] = audioread('D:\Noise filtering using Adaptive algorithm\Dataset\ref1.wav'); %
%% Algorithm's parameters
Fs = Fs_x; % Sampling rate
nwin =  512; % Window length , hamming window by default, 
noverlap = floor(0.75*nwin); % Samples of overlaps, must be smaller than nwin
nfft = 4*nwin; % Samples of FFT ,must be integer 2^n , frequncy resolution can be calculated fres=fs/nfft, 
prompt = "Filter Length: " ; % Enter filter length (must be an interger number )
N = input ( prompt ) ; %% Example:  N = 4;    
%% Algorithm steps 
%% Step 1: Compute spectral magnitudes of mixed Sound Signal and background Signal 
% 1.1: Spectrum magnitude of mixed signal, fx = spectrogram(x,nwin,noverlap,nfft) 
         
[fx , ~,~] = stft(x,Fs,'Window',hamming(nwin,'periodic'),'OverlapLength',noverlap,'FFTLength',nfft);
afx= abs(fx) ; %% spectrum magnitude
      

% 1.2: Spectrum magnitude of reference noise, fr = spectrogram(r,nwin,noverlap,nfft)     

[fr, ~,~] = stft(r,Fs,'Window',hamming(nwin,'periodic'),'OverlapLength',noverlap,'FFTLength',nfft);
afr= abs(fr) ; %% spectrum magnitude
                
 %% Step 2: Estimate channel response hh using all signal components in mixed Sound Signal and compute Spectral magnitudes of noise Signal contained in the mixed Sound Signal
[N_r, N_c]= size(afx)       ;                            %% return number of row and columns (N_C : total number of bands)
afy = zeros(N_r,N_c) ;

for i = 1:N_c
    AFR = toeplitz(afr(:,i),N) ; 
    hh = inv(transpose(AFR)*AFR)*transpose(AFR)*afx(:,i) ; 
    afy(:,i) = conv(afr(:,i),hh ); %       
end
%% Step 3: Extract silence regions in which the mixed Sound Signal Only has noise Signal components
afe = afx - afy ; %% compute residual
nafe = smooth(afe.^2) ; %% smothing afe
T = median(nafe) ;  %% applying median filter
o = nafe<T ;  %% 

%% Step 4 & Step 5 :Estimate channel response hh using only noise Signal components in mixed Sound Signal and compute Spectral magnitudes of noise Signal contained in the mixed Sound signal Step 5 & Compute spectral magnitudes of useful signal contained in the mixed Sound signal

o = reshape(o,[N_r,N_c]) ;
afy_n = zeros(N_r,N_c) ; % N_c Number of bands 
for i = 1:N_c
    afr_o=afr(:,i).*o(:,i);
    AFR = toeplitz(afr_o,N) ;
    afx_o=afx(:,i).*o(:,i) ;
    hh = inv(transpose(AFR)*AFR)*transpose(AFR)*afx_o ;  %#ok<*MINV>
    afy_n (:,i) =  conv(afr(:,i),hh ) ;            
end

afe = max(0,afx-afy_n) ;
%% Step 6: Reconstruct time domain useful signal contained in the mixed Sound signal 
S = afe.*(cos(angle(fx) +1i*sin(angle(fx)))) ; %% (phase, magnitude) to complex number 
[e,t] = istft(S,Fs,'Window',hamming(nwin,'periodic'),'OverlapLength',noverlap,'FFTLength',nfft,'Method','ola','ConjugateSymmetric',true);
%% plot, Play , write the reconstructed sound
% Time domain
figure
subplot(2,1,1)
plot(x)
title('Mixed signal')
subplot(2,1,2)
spectrogram(x,Fs)

% Frequency domain
figure
subplot(2,1,1)
plot(e)
title('Reconstructed signal')
subplot(2,1,2)
spectrogram(e,Fs)

filename = 'D:\Noise filtering using Adaptive algorithm\Dataset\Reconstructed Signals\Reconstructed_Signal.wav';
audiowrite(filename,e,Fs) ; %% save file to a given directory
soundsc(e,Fs) ; 