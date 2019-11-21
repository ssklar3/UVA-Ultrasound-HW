%  Example of use of the new Field II program running under Matlab
%
%  This example shows how a linear array B-mode system scans an image
%
%  This script assumes that the field_init procedure has been called
%
%  Example by Joergen Arendt Jensen, Version 2.0, March 22, 2011.
%  Generate the transducer apertures for send and receive
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);

f0=3e6;                  %  Transducer center frequency [Hz]
fs=100e6;                %  Sampling frequency [Hz]
c=1540;                  %  Speed of sound [m/s]
lambda=c/f0;             %  Wave length [m]
width=lambda;            %  Width of element
element_height=5/1000;   %  Height of element [m]
kerf=width/20;           %  Kerf [m]
focus=[0 0 50]/1000;     %  Fixed focal point [m]
N_elements=192;          %  Number of elements in the transducer
N_active=64;             %  Active elements in the transducer
N=10000;
%  Set the sampling frequency
x_size = 40/1000;
y_size = 10/1000;
z_size = 50/1000;
z_start = 30/1000;  %  Start of phantom surface [m];
%  Width of phantom [m]
%  Transverse width of phantom [m]
%  Height of phantom [m]
%  Create the general scatterers
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;
%  Generate the amplitudes with a Gaussian distribution
amp=randn(N,1);
%  Make the cyst and set the amplitudes to zero inside
r=5/1000;    %  Radius of cyst [m]
xc=0/1000;    %  Place of cyst [m]
zc=25/1000+z_start;
inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
amp = amp .* (1-inside);
%  Place the point scatterers in the phantom
dz=z_size/10;
for i=N-9:N
  x(i) = -15/1000;
  y(i) = 0;
  z(i) = z_start + (i-N+9)*dz;
  amp(i) = 100;
end
%  Return the variables
positions=[x y z];