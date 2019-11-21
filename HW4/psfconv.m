%% Parameters
z_foc=50e-3;             %  Range direction focal distance
fc=3e6;                  %  Transducer center frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/fc;             %  Wave length [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
num_elems=64;            %	Number of array elements
kerf=width/20;           %  Spacing between elements 
txl=kerf+width;          %  Space between signals
fieldxs=(-30e-3:txl/10:30e-3);  % X points in field
fieldzs=(30e-3:txl/10:80e-3);   % Z points in field
fs=fc/128;               % Define a sampling frequency (not critical)                                         
f=(fs:fs:10*fc);         % Define an adequate frequency range (improve resolution in time domain)
w=2*pi*f;                % Angular frequency radians
ns=length(f);            % Number of samples
bw=80;                   % Fractional bandwidth as percent
sig=bw*fc/100;           % Width of Gausian
tstep=1./max(f);         % Time steps after using Inverst FFT
t=(1:ns).*tstep;         % Define time axis
cystc=[0,50].*1e-3;      % Cyst center point x,z
cystr=5e-3;              % Cyst Radius
%% Deffinition of time delays
elsxpos=(-(txl*num_elems/2):txl:(txl*num_elems/2)); % X positions of transducer
trans_tdels=(sqrt((elsxpos.^2)+(z_foc.^2))./c);     % Time for signal from transducer to reach focus
trans_tdels=trans_tdels-max(trans_tdels);           % Built in transducer focusing delay

tdels=zeros(length(elsxpos),length(fieldxs));       % Delay matrix initialization
% calculation of time delay from each element to each x position at focus depth
for m=1:length(elsxpos)
    for i=1:length(fieldxs)
            tdels(m,i)=sqrt(((fieldxs(i)-elsxpos(m)).^2)+(z_foc.^2))./c;
    end
end
tdels=tdels-trans_tdels'; % apply built in delays to delay matrix
%% Gauss pulse claculation and summation
% Generate Gaussian pulse (frequency domain)
gauss_pulse=exp(-pi*((f-fc)/sig).^2); 
% initialize matrix of pulses
gauss_pulses=zeros(length(elsxpos),length(gauss_pulse),length(fieldxs));
% initialize matrix of pulses in time domain
gauss_t=gauss_pulses;
% Calculate pulses
for i=1:length(fieldxs)
        for m=1:length(elsxpos)
            % Calculate a pulse from each element given time delays
            gauss_pulses(m,:,i)=gauss_pulse.*exp(-j.*w.*tdels(m,i));
            % Inverse forier into time domain
            gauss_t(m,:,i)=ifft(gauss_pulses(m,:,i));  
        end
end
% Sum pulses from each transducer element at each position
gauss_t=squeeze(sum(gauss_t)); 
%% Field deffinition
% Make field of random points
pointvals=rand(length(fieldxs),length(fieldzs))-.5;
% Define cyst shadow
for i=1:length(fieldxs)
    for k=1:length(fieldzs)
        if (sqrt(((fieldzs(k)-cystc(2)).^2)+((fieldxs(i)-cystc(1)).^2)))<=cystr
            pointvals(i,k)=0;
        end
    end
end
%% Convolution of PSF with field
PSF=gauss_t(1031-20:1031+20,:); % restrict PSF to times relivant to Z range
kvs=conv2(pointvals,PSF);       % convolution
kvsr=squeeze(abs(hilbert(real(kvs))));  % Envelope the result of convolution
kvsr=kvsr./max(max(kvsr));              % Normalize
kvsr=permute(kvsr,[2 1]);               % Re-orient
%% speckle SNR calculation
A=kvsr(600:800,800:1100);   % A region of speckle
SNR=mean(A)/std(A)          % Signal to noise calculation
%% Figures
% PSF
figure
mesh(real(PSF))
xlabel('X dimension')
ylabel('Z dimension')
zlabel('Real Amplitude')
title('HomebrewPSF')
% Cyst shadow
figure
image(fieldxs*1000,fieldzs*1000,(127.*kvsr(550:1500,:)./max(max(kvsr))))
colormap(gray)
xlabel('Azimuth (mm)')
ylabel('Depth (mm)')
title('Cyst Phantom Homebew PSF')