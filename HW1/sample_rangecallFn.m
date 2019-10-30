bws=30;
fcs=10e6;
z_foc= 3e-2;
sfs=[32,64,128];
ranges=[4,8,16];
count=0;
 for n=1:length(sfs)
       for m=1:length(ranges)
           count = count+1;
           [f,att,gauss_pulse,gauss_pulse_att,gauss_t,env_gauss_t,t]=...
               GauAttFn(bws,fcs,z_foc,sfs(n),ranges(m));
            figure(1);
            subplot(3,3,count);
            plot(t,gauss_t)
            hold on
            plot(t,env_gauss_t,'--')
            title('Gaussian Pulse')
            xlabel('Time (s)')
            ylabel('Pressure - Arb Units')                                 
            title([num2str(sfs(n)),' Sampling Frequency, ', num2str(ranges(m)),...
                'X Frequency Range'])
            hold off
       end
       
 end
