%% Field call
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
%% Variables
N_elements=64;
z_foc=50e-3;             %  Range direction focal distance
fc=3e6;                  %  Transducer center frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/fc;             %  Wave length [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
num_elems=64;            %	Number of array elements
kerf=width/20;           %  Spacing between elements
txl=kerf+width;          %  Space between signals
height=5/1000;           %  hight
fieldxs=(-30e-3:txl/10:30e-3);  % X points in field
fieldzs=(30e-3:txl/10:80e-3);   % Z points in field
z=(50:txl*1000:60).*1e-3;       % Z set for PSF calc
x=(-1:txl*1000:1).*1e-3;        % X set for PSF calc
cystc=[0,50].*1e-3;      % Cyst center point x,z
cystr=5e-3;              % Cyst Radius
focus=[0,0,50e-3];       % Focus point
%% arrange points for PSF calc for field
zpointvals=zeros(length(z)*length(x),1); 
xpointvals=zpointvals;
ypointvals=zpointvals;
for i=1:length(x)
    zpointvals((length(z)*(i-1))+1:length(z)*i)=z;
    xpointvals((length(z)*(i-1))+1:length(z)*i)=x(i);
end
points=[xpointvals,ypointvals,xpointvals];
%% Field call and pressure calc
Th = xdc_linear_array (N_elements, width, height, kerf, 1, 1, focus);
[PSF,t]=calc_hp (Th,points);
PSF=PSF./max(max(PSF));
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
kvs=conv2(pointvals,PSF);               % convolution
kvsr=squeeze(abs(hilbert(kvs)));  % Envelope the result of convolution
kvsr=kvsr./max(max(kvsr));              % Normalize
kvsr=permute(kvsr,[2 1]);               % Re-orient
%% Signal to noise calculation
A=kvsr(600:800,200:400); % speckle region
SNR=mean(A)/std(A)       % SNR calculation
%% Figures
% PSF
figure
mesh(PSF)
xlabel('X dimension')
ylabel('Z dimension')
zlabel('Amplitude')
title('Field PSF')
% Cyst shadow
figure
image(fieldxs*1000,fieldzs*1000,(127.*kvsr(30:985,1:1120)./max(max(kvsr))))
colormap(gray)
xlabel('Azimuth (mm)')
ylabel('Depth (mm)')
title('Cyst Phantom Field PSF')