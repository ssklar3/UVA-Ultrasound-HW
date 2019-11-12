%% Field call
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
%% Variables
a=5e-3;                 % Half x width
lam=1e-3;               % Wavelength
zend=((a^2)./lam)*1e3;        % End Z point
z=(.1:.1:zend).*1e-3;   % Z field
y=0;                    % Y field
x=(0:.1:1.5).*5e-3;     % X field
%% Field Parameters
frq=1e10;
set_sampling(frq);
center=[0,0,0];
focus=[0,0,12.5e-3];
tmax=4e-5;
c=1540;
%% Field call and pressure calc
Th = xdc_piston (5e-3, 5e-6);  % Transducer definition

p=zeros(length(x),length(z));            % Pressure array initialization
for i=1:length(z)
    i
    for n= 1:length(x)
        [h,t]=calc_h (Th,[x(n),y,z(i)]); % Transfer function calc
        t=[t:1/frq:tmax];                % Define times for sin wave
        wav=sin(2*pi*t*c./lam);          % creation of sin wave
        hwc=conv(h,wav);                 % convolution
        hwc=diff(hwc);                   % differentiation
        prS=floor((length(hwc)./2)-(length(hwc)./5));   % Lower bounding for amplitude check
        prE=ceil((length(hwc)./2)+(length(hwc)./5));    % Higher bounding for amplitude check
        p(n,i)=max(hwc(prS:prE))-min(hwc(prS:prE));          % Calculate pressure
    end
end
p=p./max(max(p)); % Pressure Normalization
%% plot
figure
    mesh(z/((a^2)/lam),x/(5e-3),p);
    xlabel('Depth (a^2/\lambda)')
    ylabel('Azimuth (a)')
    zlabel('Normalized Pressure')
    title('Circular Mesh for a/\lambda=5')