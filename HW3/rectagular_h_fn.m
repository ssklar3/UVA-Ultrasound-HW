function [p,h,t,hwc] = rectagular_h_fn (x,y,z,lam,rect)
path(path,'/Users/samuelsklar/Documents/Classwork/Ultrasound/Field_II_ver_3_24_mac');
field_init(-1);
frq=2e11;
set_sampling(frq);
a=5e-3;
b=7.5e-3;
center=[0,0,0];
focus=[0,0,50e-3];
tmax=4e-5;
c=1540;
A=5e-3;
B=7.5e-3;
% pitch=.1e-3;
% xpts=-2*a:pitch:2*a;
% ypts=-2*b:pitch:2*b;
% noels=(length(xpts)-1)*(length(ypts)-1);
% rect=zeros(noels,19);
% rect(:,1)=1;
% rect(:,14)=1;
% rect(:,15)=pitch;
% rect(:,16)=pitch;
% count=0;
% for j=1:length(ypts)-1
%     for k=1:length(xpts)-1
%         count=count+1;
%         rect(count,2:3)=[xpts(k),ypts(j)];
%         rect(count,5:6)=[xpts(k+1),ypts(j)];
%         rect(count,8:9)=[xpts(k+1),ypts(j+1)];
%         rect(count,11:12)=[xpts(k),ypts(j+1)];
%         rect(count,17:18)=[mean(xpts(k:k+1)),mean(ypts(j:j+1))];
%     end
% end
% 
% % rect=[1,-a,-b,0, a,-b,0, a,b,0, -a,b,0, 1,2*a,2*b,0,0,0,];
% 
Th = xdc_rectangles (rect, center, focus);
[h,t]=calc_h (Th,[x,y,z]);
tend=t+(length(h)/frq);
%h=[(0:1/frq:t)*0,h',(tend:1/frq/tmax)];
% plot(h)
t=[t:1/frq:tmax];
% dh=diff(h);                % differentiation of transfer function
% wav=sin(2*pi*t*c./lam);    % creation of sin wave
% hwc=conv(dh,wav);               % convolution


wav=sin(2*pi*t*c./lam);    % creation of sin wave
hwc=conv(h,wav);               % convolution
hwc=diff(hwc);
prS=floor((length(hwc)./2)-(length(hwc)./5));   % Lower bounding for amplitude check
prE=ceil((length(hwc)./2)+(length(hwc)./5));    % Higher bounding for amplitude check
p=max(hwc(prS:prE))-min(hwc(prS:prE));          % Calculate pressure\

end