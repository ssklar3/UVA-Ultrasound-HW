%% Circular piston paramaters
c=1540;                 % Speed of sound in body (m/s)
z=25e-3;                % Range (m)
a=5e-3;                 % Piston Radius (m)
frq=1e11;               % Freqency Hz
r=[6,5,4.5,2.5,0.5,0]*1e-3; % Radial distance of point from axis
ts=cell(length(r),1);   % Bin for times within the detection time frame
tcon=ts;                % Bin for times within the h=1 time frame
hs=ts;                  % Bin for transfer fn values
tpre=ts;                % Bin for times before first detection
tpost=ts;               % Bin for times after detection frame
tmax=1.6*(sqrt((z^2)+((a+max(r)).^2))./c);  % End of transfer fn time window
t=zeros(length(r),floor(tmax*frq));         % Finnal time array
h=t;                                        % Finnal transfer functino array
figure
for i=1:length(r)       % Loop through positions
    tff=0;
    %% h calculation
    if a>r(i)           % Define time of first detection if point within piston radius
         tff=z./c;
         tf=sqrt((z.^2)+(r(i)-a).^2)./c;
    else                % Define time of first detection if point out of piston radius
        tf=sqrt((z.^2)+(r(i)-a).^2)./c;
    end
    te=sqrt((z^2)+((a+r(i)).^2))./c; % Define time of last detection 
    ts{i}=[tf:(1/frq):te];           % Array of times durring detection
    ts{i}=ts{i}(2:end-1);              % Remove overlap between time cells
    % h calculating line 
      hs{i}=(1/pi).*acos(((((c.^2).*(ts{i}.^2))-(z.^2)+(r(i).^2)-(a.^2))...
        ./(2.*r(i).*((c.^2.*ts{i}.^2)-z.^2).^(.5))));
    %% array construction 
    if a>r(i)   
        tcon{i}=[tff:(1/frq):tf];       % Array of times where detection only in apature
        tpre{i}=[0:(1/frq):tff];        % Array of times before detection
        tpre{i}=tpre{i}(1:end-1);
    else
        tpre{i}=[0:(1/frq):tf];          % Array of times before detection
        tcon{i}=[];
    end
    
    tpost{i}=[te:(1/frq):tmax];      % Array of times after detection
    
    % Checks to see if the right number of times exist if not adds another
    % point on the end or removes one
    % (I could get it to work wihtout this but it is more robust this way)

    while length(tpre{i})+length(tcon{i})+length(ts{i})+length(tpost{i}) < length(t(i,:))
        tpost{i}=[tpost{i},tpost{i}(end)+(1/frq)];
    end
    while length(tpre{i})+length(tcon{i})+length(ts{i})+length(tpost{i}) > length(t(i,:))
        tpost{i}=tpost{i}(1:end-1);
    end
    
    t(i,:)=[tpre{i},tcon{i},ts{i},tpost{i}]; % Combine time cells into finnal array  
%% output of transfer function
    h(i,length(tpre{i})+length(tcon{i})+1:length(tpre{i})+length(tcon{i})+length(ts{i}))=hs{i};
% Defines points where h=1
    h(i,length(tpre{i})+1:length(tpre{i})+length(tcon{i}))=1;
    %% plot
    subplot(length(r),1,i)
    plot(t(i,:),h(i,:))
    xlabel('Time (s)')
    ylabel('h(x,t)/c')
    axis([1.6e-5 1.8e-5 0 1.2])
    title(['r/a = ',num2str(r(i)/a)])
end