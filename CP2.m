clear all; close all; clc;
load 'CP2_SoundClip.mat'

Fs = 44100; %sample rate of the sound clip

S = y'; % transposes in order to have consistent dimensions
w = length(y)/4; % break the spectogram up into four time windows, otherwise it
% will be too big for MATLAB/autograder to run.

S1 = S((1-1)*w+1:1*w); % this will isolate the correct window
S2 = S((2-1)*w+1:2*w); % this will isolate the correct window
S3 = S((3-1)*w+1:3*w); % this will isolate the correct window
S4 = S((4-1)*w+1:4*w); % this will isolate the correct window

L = length(S1)/Fs; % length of each window in seconds
n = length(S1);  % number of elements in each window
t = [0:1/Fs:L - 1/Fs]; % t in sec. relative to the start of the window
tau = 0:0.1:L; % discretization for the Gabor transform
k = 2*pi*(1/L/2)*[0:n/2-1 -n/2:-1]; % discretization in frequency space
ks = fftshift(k); % gotta shift them freqs.

Sgt_spec = zeros(length(ks),length(tau)); % initializing the function for the spectrogram

%Gabor Transform Parameters
a = 400; %this will give you the correct width, so exp(-a(...))
range = [1:1800]; %use this when finding the max and index of the transformed Gabor filtered signal
% i.e., max(TransformedSignal_GaborFiltered(range))
%%

% Repeat this part for S1, S2, S3, and S4.
% For this part we are going to make a spectrogram, but only for the freqs.
% of interest.  So first we'll do a Gabor transform, then we'll filter
% around our peak freq. in the regime of interest, and then we'll look at
% the spectrogram of that function.

% Gabor transform each S, just like we did in the lecture
% Week4_Spectrograms.m lines 141 to 146.
% You'll have to add code between line 144 and Sgt_spec at line 145.



% After line 144 find the index of your peak frequency (in absolute value)
% within the range of interest (this range is very forgiving so you don't
% have to match the autograder exactly, just use your judgement based on
% the figure in the assignment statement).  Then make a filter centered
% about the peak frequency.  Filter the Gabor transformed function.
% this is the function you will use in line 145 from the lecture to find
% your Sgt_spec.

% apply Gabor Transform to achieve spectogram, only keeping the peak
% frequencies 
%for j = 1:length(tau)
%   g = exp(-a*(t - tau(j)).^2); % Window function
%   % apply gauss filter to reading (time amp space)
%   Sg = g.*S1;
%   %apply FFT (k, amp space)
%   Sgt = fft(Sg);
%
%   % find index of peak frequency
%   [maxVal,maxIndex] = max(Sgt(range));
%   
%   % create gauss filter in frequency space centered around peak frequency
%   g2 = exp((-1/L)*(abs(k)-k(maxIndex)).^2);
%
%   % apply filter to FFT
%   Sgt_filt2 = g2.*Sgt;
%
%   Sgt_spec(:,j) = fftshift(abs(Sgtip_filt2)); % We don't want to scale it
%end
%
%A1 = Sgt_spec;


for j = 1:length(tau)
   g = exp(-a*(t - tau(j)).^2); % Window function
   % apply gauss filter to reading (time amp space)
   Sg = g.*S2;
   %apply FFT (k, amp space)
   Sgt = fft(Sg);

   % find index of peak frequency
   [maxVal,maxIndex] = max(Sgt(range));
   
   % create gauss filter in frequency space centered around peak frequency
   g2 = exp((-1/L)*(abs(k)-k(maxIndex)).^2);

   % apply filter to FFT
   Sgt_filt2 = g2.*Sgt;

   Sgt_spec(:,j) = fftshift(abs(Sgt)); % We don't want to scale it
end

A2 = Sgt_spec;

%for j = 1:length(tau)
%   g = exp(-a*(t - tau(j)).^2); % Window function
%   % apply gauss filter to reading (time amp space)
%   Sg = g.*S3;
%   %apply FFT (k, amp space)
%   Sgt = fft(Sg);
%
%   % find index of peak frequency
%   [maxVal,maxIndex] = max(Sgt(range));
%   
%   % create gauss filter in frequency space centered around peak frequency
%   g2 = exp((-1/L)*(abs(k)-k(maxIndex)).^2);
%
%   % apply filter to FFT
%   Sgt_filt2 = g2.*Sgt;
%
%   Sgt_spec(:,j) = fftshift(abs(Sgt_filt2)); % We don't want to scale it
%end
%
%A3 = Sgt_spec;
%for j = 1:length(tau)
%   g = exp(-a*(t - tau(j)).^2); % Window function
%   % apply gauss filter to reading (time amp space)
%   Sg = g.*S4;
%   %apply FFT (k, amp space)
%   Sgt = fft(Sg);
%
%   % find index of peak frequency
%   [maxVal,maxIndex] = max(Sgt(range));
%   
%   % create gauss filter in frequency space centered around peak frequency
%   g2 = exp((-1/L)*(abs(k)-k(maxIndex)).^2);
%
%   % apply filter to FFT
%   Sgt_filt2 = g2.*Sgt;
%
%   Sgt_spec(:,j) = fftshift(abs(Sgt_filt2)); % We don't want to scale it
%end
%
%A4 = Sgt_spec;




% Save Sgt_spec as variable A1 after your for loop.  Repeat for S2, S3, and
% S4 and save those Sgt_spec as A2, A3, and A4.  You don't have to rewrite
% the code, just copy and paste and use the respective S's, or write a for
% loop that iterates through S1 to S4.

%A1 =     % Shape:  484560x110 double
%A2 =     % Shape:  484560x110 double
%A3 =     % Shape:  484560x110 double
%A4 =     % Shape:  484560x110 double


% Plot of spectrogram for each window (for the report, not for autograder)
% just like we did in the lecture, but change your ylim to be in the range
% of interest for the sound clip.

%figure(1)
%pcolor(tau,ks,A1)
%shading interp
%set(gca,'ylim',[0 500],'Fontsize',16)
%colormap(hot)
%colorbar
%xlabel('time (t)'), ylabel('frequency (k)')
%title('S1 Spectogram')
%
figure(2)
hold on;
pcolor(tau,ks,A2)
shading interp
set(gca,'ylim',[0 500],'Fontsize',16)
colormap(hot)
colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title('S2 Spectogram')


%figure(3)
%hold on;
%pcolor(tau,ks,A3)
%shading interp
%set(gca,'ylim',[0 500],'Fontsize',16)
%colormap(hot)
%colorbar
%xlabel('time (t)'), ylabel('frequency (k)')
%title('S3 Spectogram')


%figure(4)
%hold on;
%pcolor(tau,ks,A4)
%shading interp
%set(gca,'ylim',[0 500],'Fontsize',16)
%colormap(hot)
%colorbar
%xlabel('time (t)'), ylabel('frequency (k)')
%title('S4 Spectogram')




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Isolate the bassline
S = y'; % transposes in order to have consistent dimensions
L = length(S)/Fs; % total length in sec. 
n = length(S); % total number of elements in S
t = [0:1/Fs:L - 1/Fs]; % time discretization
k = (1/L)*[0:n/2-1 -n/2:-1]; % freq. discretization
kt = 2*pi*(1/L/2)*[0:n/2-1 -n/2:-1]; % freq discretization 
ks = fftshift(kt); % shift freq discretization so from least to greatest 

% Take the Fourier transform of S, and in freq. space isolate all freqs. 
% (in absolute value) that you determine should be part of the baseline 
% according to spectrogram (or also just by listening); that is, all points
% in the transformed function not within the frequency range you determined
% should be set to zero (kind of like a Shannon filter, but simpler than
% what we did in lecture).
% You may have to do this part a few times with different thresholds to get
% it right.

% take fourier transf of S
St = fftshift(fft(S));
% find max value for plot
maxF = max(St);
figure(2)
plot(ks,abs(St)/maxF);

%frequencies that should be part of bass around 275 Hz
% filter frequencies, set all values for frequenices not in 250-300 to be 0
% find index of frequenicies (0 is ocated at size/2+1)
f1_index = size(ks,2)/2+1+250;
f2_index = size(ks,2)/2+1+300;

% create new matrix of all zeroes
Sfilt = zeros(size(St));
% copy over values in found range in St to new matrix 
Sfilt(1,f1_index:f2_index) = St(1,f1_index:f2_index);


figure(3)
plot(ks,abs(Sfilt))

%figure(2)
%plot(ks,St)
%shading interp
%set(gca,'ylim',[0 500],'Fontsize',16)
%colormap(hot)
%colorbar
%xlabel('time (t)'), ylabel('frequency (k)')


% After thresholding the transformed function, take the inverse transform
% of the thresholded function and save it as A5.

A5 = ifft(Sfilt)';    %Shape:  1938240x1 double

%Play sound (not for autograder)


%Plot the amplitude S over time (for the report, not for the autograder)
figure(4)
plot(t,A5);
xlabel('time (t)'), ylabel('Amplitude')
title('Isolated Bass Sound Clip')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Isolate the guitar

%Same exact process as the baseline above, but you'll have to be more
%careful about the frequency range.
S = y'; %reinitialize the S from the previous part above.

% take fourier transf of S
St = fftshift(fft(S));

%frequencies that should be part of bass around 275 Hz
% filter frequencies, set all values for frequenices not in 260-480 to be 0
% find index of frequenicies (0 is ocated at size/2+1)
f1_index = size(ks,2)/2+1+460;
f2_index = size(ks,2)/2+1+480;

% create new matrix of all zeroes
Sfilt = zeros(size(St));
% copy over values in found range in St to new matrix 
Sfilt(1,f1_index:f2_index) = St(1,f1_index:f2_index);

A6 =  ifft(Sfilt)';    %Shape:  1938240x1 double

%Play sound (not for autograder)


%Plot the amplitude S over time (for the report, not for the autograder)
figure(5)
plot(t,A6);
xlabel('time (t)'), ylabel('Amplitude')
title(['Isolated Guitar Sound Clip'])
