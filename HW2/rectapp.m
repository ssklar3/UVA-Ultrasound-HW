function [h,t,p]=rectapp(x,y,z,lam,tmax)
%% rectangular paramaters
c=1540;                 % Speed of sound in body (m/s)
a=5e-3;                 % short axis (m)
b=7.5e-3;               % long axis (m)
xset= zeros(4,1);       % Set of box x axis lengths
yset= zeros(4,1);       % set of box y axis lengths
sset= zeros(4,1);       % Set of box short axis lengths
lset= zeros(4,1);       % set of box long axis lengths
frq=1e10;               % Freqency Hz
critTs=zeros(4);        % Critical time bin
tpre=cell(4,1);   % Bin for times before the detection time frame
toneset=tpre;     % Bin for when h=1
tshort=tpre;      % Bin for short detection before long
tlong=tpre;       % Bin between long and end
tpost=tpre;       % Bin for times after detection time frame
hs=cell(4,5);                  % Bin for transfer fn values
t=zeros(1,floor(tmax*frq));    % Finnal time array
h=t;                           % Finnal transfer functino array
hsq=zeros(4,floor(tmax*frq));  % H components from each rectangle
sign=ones(4,1);
%% Class of position
if abs(x)<=a
    if abs(y)<=b
        class=1;        % Inside apature
    else
        class=2;        % outside in y dimension
    end
else
    if abs(y)<=b
        class=3;        % outside in x dimension
    else
        class=4;        % outside in both dimensions
    end
end
%% Loop through rectangles
for i=1:4
    %% Define regtangle dimensions and sign of opperation
    if class==1
        xset(i)=a-abs(x).*(-1)^i;
        yset(i)=b-abs(y).*(-1)^ceil(i/2);
    end
    if class==2
        xset(i)=a-abs(x).*(-1)^i;
        if i<=2
            yset(i)=abs(y)-b;
            sign(i)=-1;
        else
            yset(i)=2*b;
        end
    end
    if class==3
        yset(i)=b-abs(y).*(-1)^i;
        if i<=2
            xset(i)=abs(x)-a;
            sign(i)=-1;
        else
            xset(i)=2*a;
        end
    end
    if class==4
        if i==1
            yset(i)=2*b;
            xset(i)=2*a;
        end
        if i==2
            yset(i)=abs(y)-b;
            xset(i)=2*a;
            sign(i)=-1;
        end
        if i==3
            xset(i)=abs(x)-a;
            yset(i)=2*b;
            sign(i)=-1;
        end
        if i==4
            yset(i)=abs(y)-b;
            xset(i)=abs(x)-a;
        end
    end
    %% Order sets by size
    if xset(i)>yset(i)
        sset(i)=yset(i);
        lset(i)=xset(i);
    else
        sset(i)=xset(i);
        lset(i)=yset(i);        
    end
    %% Time set up
    critTs(i,1)=z/c; % Define first critical time
    critTs(i,2)=sqrt((z^2)+(sset(i)^2))/c; % Define second critical time
    critTs(i,3)=sqrt((z^2)+(lset(i)^2))/c; % Define third critical time
    critTs(i,4)=sqrt((z^2)+(sset(i)^2)+(lset(i)^2))/c; % Define fourth critical time
    tpre{i}=[0:(1/frq):critTs(i,1)]; % Times before detection
    toneset{i}=[critTs(i,1):(1/frq):critTs(i,2)]; % Times where responce is one
    toneset{i}=toneset{i}(2:end); % Shorten front to avoid overlap
    tshort{i}=[critTs(i,2):(1/frq):critTs(i,3)]; % Times after short edge is reached
    tshort{i}=tshort{i}(2:end); % Shorten front to avoid overlap
    tlong{i}=[critTs(i,3):(1/frq):critTs(i,4)];% Times after long edge is reached
    tlong{i}=tlong{i}(2:end-1); % Shorten front and end to avoid overlap
    tpost{i}=[critTs(i,4):(1/frq):tmax]; % Times after detection
    %% make sure there are no size artifacts
    while length(tpre{i})+length(toneset{i})+length(tshort{i})...
        +length(tlong{i})+length(tpost{i})>length(t)
        tpost{i}=tpost{i}(1:end-1);
    end
    while length(tpre{i})+length(toneset{i})+length(tshort{i})...
        +length(tlong{i})+length(tpost{i})<length(t)
        tpost{i}=[tpost{i},0];
    end
    %% set of h calculations for each time frame
    hs{i,1}=tpre{i}*0;
    hs{i,2}=ones(1,length(toneset{i})).*(pi/2);
    hs{i,3}=(pi/2)-acos(sset(i)./sqrt((c^2).*(tshort{i}.^2)-(z^2)));
    hs{i,4}=(pi/2)-acos(sset(i)./sqrt((c^2).*(tlong{i}.^2)-(z^2)))...
        -acos(lset(i)./sqrt((c^2).*(tlong{i}.^2)-(z^2)));
    hs{i,5}=tpost{i}*0;

    hsq(i,:)=[hs{i,1},hs{i,2},hs{i,3},hs{i,4},hs{i,5}];

    h=h+(hsq(i,:).*sign(i));
end

t=[tpre{i},toneset{i},tshort{i},tlong{i},tpost{i}];
h=h./(2.*pi);
%% pressure calc
dh=diff(h);                % differentiation of transfer function
wav=sin(2*pi*t*c./lam);    % creation of sin wave
hwc=conv(dh,wav);               % convolution
prS=floor((length(hwc)./2)-(length(hwc)./5));   % Lower bounding for amplitude check
prE=ceil((length(hwc)./2)+(length(hwc)./5));    % Higher bounding for amplitude check
p=max(hwc(prS:prE))-min(hwc(prS:prE));          % Calculate pressure
end