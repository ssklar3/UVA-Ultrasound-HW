bws=30;
fcs=[5e6,10e6];
z_foc= 50e-2;
sfs=64;
ranges=8;
rez=zeros(length(fcs),1);
 for n=1:length(fcs)
           [f,att,gauss_pulse,gauss_pulse_att,gauss_t,env_gauss_t,t]=...
               GauAttFn(bws,fcs(n),z_foc,sfs,ranges);
            st=spline(env_gauss_t(10:40*n),t(10:40*n),.1);
            ed=spline(env_gauss_t(50*n:70*n),t(50*n:70*n),.1);
            rez(n) = (ed-st)*1540.*2
 end