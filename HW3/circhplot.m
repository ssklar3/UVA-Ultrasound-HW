%% Field call
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
%% Variables
a=5e-3;                 % Half x width
x=[1.2,1,0.9,0.5,0.1,0].*5e-3;    % X points
y=zeros(1,length(x));    % Y points
z=(50).*1e-3;           % Z point
z=z*ones(1,length(x));
%% Field Parameters
frq=1e11;
set_sampling(frq);
center=[0,0,0];
focus=[0,0,50e-3];
tmax=4e-5;
c=1540;
%% Field call
Th = xdc_piston (5e-3, 5e-6);  % Transducer definition
[h,t]=calc_h (Th,[x',y',z']); % Transfer function calc
h=h./max(max(h));
%% Add padding on front end
hfront=zeros(ceil(((t-(3e-5))*frq)),6);
tf=(3e-5:1/frq:t+length(h)/frq);
h=[hfront;h];
%% plot
figure
for i=1:length(x)
    subplot(length(x),1,i)
    plot(tf,h(:,i));
    xlabel('Time (s)')
    ylabel('h(x,t)/c')
    axis([3.24e-5 3.35e-5 0 1.2])
    title(['x/a = ',num2str(x(i)/5e-3)])
end