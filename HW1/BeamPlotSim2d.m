z_foc=50e-3;                                % Range direction focal distance
%theta_steer_vert_deg=0;
%theta_steer=(90-theta_steer_vert_deg)*pi/180;                       % Steer angle - will use this later
theta_steer=(0)*pi/180;
num_elems=32;                              % Number of array elements
fc=10e6;                                    % Center frequency
pitch=(.075e-3);                         
vel=1540;                                   % Speed of sound - all units MKS
fs=fc/64;                                   % Define a sampling frequency (not critical)                                         
f=[fs:fs:8*fc];                             % Define an adequate frequency range (improve resolution in time domain)
i_cen=round(fc/fs);                         % Index of center frequency for single frequency calc
w=2*pi*f;                                   % Angular frequency radians
ns=length(f);                               % Number of samples
tdel=0e-6;                                  % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
tdel2=1.0e-6;                               % Use a fixed time offset so that initial waveform is 0 for t<0 (not critical)
bw=30;                                      % Fractional bandwidth as percent
bw2=80;                                     % Fractional bandwidth as percent
sig=bw*fc/100;                              % Width of Gausian
gauss_pulse=exp(-pi*((f-fc)/sig).^2);       % Generate Gaussian pulse (frequency domain)
gauss_pulse=gauss_pulse.*exp(-j*w*tdel);    % Apply time delay so 0 for t<0
%att= .00005*f;
att=0;
%%
weight=ones(num_elems,1);                  % Define the weighting function (you can change this)
%weight=hann(num_elems);

x_pts=[0:0.2:5].*1e-3;
y_pts=[0:0.2:5].*1e-3;
z_pts=[50].*1e-3;
%num_elems=64;
x_elem=zeros(num_elems,1);
y_elem=zeros(num_elems,1);
foc_del=zeros(num_elems);
for i=1:num_elems
    x_elem(i)=((i-1)-(num_elems-1)./2).*pitch;  % Calculate locations of each array element]
    for k=1:num_elems
        y_elem(k)=((k-1)-(num_elems-1)./2).*pitch;
        foc_del(i,k) =  (sqrt(((x_elem(i))^2)+(y_elem(k)^2)+(z_foc^2)))./vel;
    end
end

field_val= zeros(length(x_pts),length(y_pts));
field_val_fc= zeros(length(x_pts),length(y_pts));
field_db= zeros(length(x_pts),length(y_pts));
field_fc_db= zeros(length(x_pts),length(y_pts));

    for i=1:length(x_pts)
        for n=1:length(y_pts)
        sum_pulse=zeros(size(gauss_pulse));                 % Initialize sum of waveforms to zero before starting loop
            for k=1:num_elems
                for m=1:num_elems
                r=sqrt(((x_elem(k)-x_pts(i))^2)+((y_elem(m)-y_pts(n))^2)+(z_pts^2));
                prop_del = (r./vel);
                  sum_pulse=sum_pulse+weight(k).*10.^(-att.*r./20).*gauss_pulse.*... % Sum pulse taking account of current pulse, propagation delay - focusing delay
                    exp(-sqrt(-1)*w.*(prop_del-foc_del(k,m)));
                end
            end
        
        sum_pulse_t=real(ifft(sum_pulse));                  % Convert to time domain

        field_val(i,n)= max(abs(hilbert(sum_pulse_t)));     % Field values in beamplot - Hilbert finds envelope of sinsoidal fn
        end
    end
     field_val(:,:)=field_val(:,:)./max(max(field_val(:,:)));     % Normalize by peak value at each range
     field_db(:,:)=20*log10(field_val(:,:));                 % Convert to dB

%%
field_db=[flip(field_db(:,1:end-1),2),field_db];
field_db=[flip(field_db(1:end-1,:),1);field_db];
x_ptsp=[-5:0.2:5].*1e-3;
y_ptsp=[-5:0.2:5].*1e-3;
figure;
contour(y_ptsp*1e3,x_ptsp*1e3,field_db,[-30 -20 -15 -12 -9 -6 -3]);
title('Contour Field')
xlabel('Depth (mm)')
ylabel('Azimuth (mm)')
figure
mesh(y_ptsp*1e3,x_ptsp*1e3,field_db);
title('Mesh Field')
xlabel('Depth (mm)')
ylabel('Azimuth (mm)')
zlabel('Signal (dB)')