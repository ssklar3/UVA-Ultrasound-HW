function SSE = SSEleptinuptake(p)
% Data
% Plasma leptin concentration
S = [0.75; 1.18; 1.47; 1.61; 1.64; 5.26; 5.88; 6.25; 8.33; 10.0; 11.11; ...
20.0; 21.74; 25.0; 27.77; 35.71];
% Renal Leptin Uptake
R = [0.11; 0.204; 0.22; 0.143; 0.35; 0.48; 0.37; 0.48; 0.83; 1.25; ...
0.56; 3.33; 2.5; 2.0; 1.81; 1.67];
% model variables Rmax = p(1);
Km = p(2);
Rmax =p(1);
SSE = sum((R - Rmax*S./(Km + S)).^2);
end