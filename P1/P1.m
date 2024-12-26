clear all
close all
clc

R0 = 1;
Kp = 25/2;
Kd = 1;
Gs = 1;
Ga = 0.095;
Gr = 1;
Gd = 1;
Gf = 1;
s = tf('s');
nu = 0;
p = 1;

% nominal model definition
%PUI
Ku=[9;16];           Kn=(Ku(1)+Ku(2))/2; 
p1u=[0.55;1.05];     p1n=(p1u(1)+p1u(2))/2;
p2u=[1.9;3.1];       p2n=(p2u(1)+p2u(2))/2; 

%nominal model for the plant
Gpn = minreal(zpk(Kn/(s*(1+s/p1n)*(1+s/p2n))))

%S2
rho_r = 1.5e-1;
S_star0_S2 = rho_r / abs(Kd*R0);

%S3
Da0 = 5.5e-3;
rho_a = 1.5e-2;
S_star0_S3 = rho_a / abs(Kp*Da0);

%S4
rho_p = 5e-4;
ap = 2e-2;
MS_LF = rho_p/ap;
MS_LF_dB = 20*log10(MS_LF);
wp_plus = 0.02;

%S5
rho_s = 5e-4;
as = 1e-1;
MT_HF = (rho_s*Gs) / as;
MT_HF_dB = 20*log10(MT_HF);
ws_minus = 40;

%S6-S8
zeta_II = abs(log(0.1)) / sqrt(pi^2+(log(0.1))^2);
Sp0 = 2*zeta_II*sqrt( 2+4*zeta_II^2+2*sqrt(1+8*zeta_II^2) ) / ( sqrt(1+8*zeta_II^2)+4*zeta_II^2-1 );
Sp0_dB = 20*log10(Sp0);
Tp0 = 1 / ( 2*zeta_II*sqrt(1-zeta_II^2) );
Tp0_dB = 20*log10(Tp0);
wc1 = (1/3) * (1/(sqrt(1-zeta_II^2))) * (pi-acos(zeta_II)) * sqrt(sqrt(1+4*zeta_II^4)-2*zeta_II^2);
wc2 = (1/12) * (-log(0.05)/zeta_II) * sqrt(sqrt(1+4*zeta_II^4)-2*zeta_II^2);
wn_II1 = (1/(3*sqrt(1-zeta_II^2))) * (pi-acos(zeta_II));
wn_II2 = -log(0.05) / (12*zeta_II);
wc = max(wc1, wc2);
wn_II = max(wn_II1, wn_II2);

% % ------------------Weighting function on S------------------------------
S_star0 = 0.15;
S_to0=S_star0*s^(nu+p); 
S_II=(s*(s/wn_II^2+2*zeta_II/wn_II))/(1+(2*zeta_II/wn_II)*s+s^2/wn_II^2);
yline(Sp0_dB), yline(MS_LF_dB), xline(wp_plus)

WS_inv=0; 
hold on, grid on, zoom on
bodemag(S_to0,'r--', S_II, 'b--');
a = S_star0
z1 = 0.001
p1 = 0.0086
w1 = 0.78

zeta_WS = 0.81;
z2 = (a*p1*w1^2)/(z1*Sp0)
WS_inv = ( a*s*(1+s/z1)*(1+s/z2) )/( (1+s/p1)*(1+(2*zeta_WS/w1)*s+s^2/w1^2) );
WS_inv = zpk(WS_inv);
bodemag(WS_inv);

% % -----------------Weighting function on T---------------------------
figure(2); 
hold on, grid on, zoom on
T_II=1/(1+(2*zeta_II/wn_II)*s+s^2/wn_II^2);

yline(Tp0_dB), xline(ws_minus), yline(MT_HF_dB)
bodemag(T_II, 'b--')

wT = ws_minus * 10^((MT_HF_dB-Tp0_dB)/40)
WT_inv = Tp0/(1+(1.414/wT)*s+s^2/wT^2);
WT_inv = zpk(WT_inv);

bodemag(WT_inv) 

% % -----------------Design of Wu--------------------------------
 nsamp = 10;
 num = 500;
 omega = logspace(-3,3, num);
% gridK = linspace(Ku(1),Ku(2),nsamp);
% gridp1 = linspace(p1u(1), p1u(2),nsamp); 
% gridp2 = linspace(p2u(1),p2u(2),nsamp); 
% matDelta = zeros(nsamp^3, num); cnt = 1;
% 
% % cloud of transfer function
% figure("Name", "Weighting function Wu") 
% 
% for ii = 1:nsamp
%     for jj = 1:nsamp
%         for kk = 1:nsamp
%             K_ = gridK(ii);
%             p1_ = gridp1(jj);
%             p2_ = gridp2(kk);
% 
%             Gp = minreal(zpk(K_ / (s*(1+s/p1_)*(1+s/p2_))));
%             Delta = (Gp/Gpn)-1;
%             [magDelta, ~] = bode(Delta, omega);
%             magDelta = squeeze(magDelta);
%             matDelta(cnt, :) = magDelta';
%             cnt = cnt+1;
%             semilogx(omega, 20*log10(magDelta));
%             hold on, grid on;
%         end
%     end
% end
% 
% maxDelta = max(matDelta)';
% Wu_mag = vpck(maxDelta, omega);
% fit = fitmag(Wu_mag);
% [A, B, C, D] = unpck(fit);
% [numWu, denWu] = ss2tf(A, B, C, D);
% Wu = tf(numWu, denWu);
% 
% figure("Name", "Weighting function Wu") 
% for ii = 1:nsamp
%     for jj = 1:nsamp
%         for kk = 1:nsamp
%             K_ = gridK(ii);
%             p1_ = gridp1(jj);
%             p2_ = gridp2(kk);
% 
%             Gp = minreal(zpk(K_ / (s*(1+s/p1_)*(1+s/p2_))));
%             Delta = (Gp/Gpn)-1;
%             [magDelta, ~] = bode(Delta, omega);
%             magDelta = squeeze(magDelta);
%             matDelta(cnt, :) = magDelta';
%             cnt = cnt+1;
%             semilogx(omega, 20*log10(magDelta));
%             hold on, grid on;
%         end
%     end
% end
% 
% [magWu,~]=bode(Wu,omega); 
% magWu=squeeze(magWu);
% semilogx(omega,20*log10(magWu),'LineWidth',2, 'Color','red');

Wu= minreal(zpk((1.083*s^3 + 3.543*s^2 + 2.283*s + 0.3379)/...
    (s^3 + 4.897*s^2 + 5.653*s + 1.208)));

%% LMI H_infinity design

% modification of W1
wc = 1;
lambda = 0.01*wc;
W1 = minreal(zpk(inv(WS_inv)), 1e-3); % WS
W1mod = minreal(zpk( W1*(s/(s+lambda))^(nu+p) ), 1e-3);

% choice of W2
% hold off
% figure("Name", "Wu vs WT")
% [magWu, ~] = bode(Wu, omega);
% magWu = squeeze(magWu);
% semilogx(omega, 20*log10(magWu), 'LineWidth', 1, 'Color', 'red');
% hold on, grid on, zoom on
% WT = minreal(zpk(inv(WT_inv)), 1e-3);
% [magWT, ~] = bode(WT, omega);
% magWT = squeeze(magWT);
% semilogx(omega, 20*log10(magWT), 'LineWidth', 1, 'Color', 'blue');

% W2 = WT
W2mod = tf(1, Tp0);

% state space representation of the generalized plant for RS and NP
[Am,Bm,Cm,Dm] = linmod('P1_GeneralizedPlant_RS_NP');
M = ltisys(Am,Bm,Cm,Dm);

%I insert the zeros I have removed
M = sderiv(M,2,[1/wT 1]);
M = sderiv(M,2,[1/wT 1]);
[gopt, Gcmod] = hinflmi(M, [1 1], 0, 0.01, [0 0 0]);
[Ac, Bc, Cc, Dc] = ltiss(Gcmod);
Gcmod = ss(Ac, Bc, Cc, Dc);
Gcmod = minreal(zpk(Gcmod),1e-4)
[zGc, pGc, kGc] = zpkdata(Gcmod,'v');
Gc = minreal(zpk( Gcmod*((1+s/5453)*(s+0.01))/(s+5.304e-6) ), 1e-3)

figure("Name", "Nichols diagram");
myngridst(Tp0,Sp0); 
hold on,grid on

%% Check of performances

%check of nominal stability  
Ln = minreal(zpk(Gc*Gpn*Gs*Ga*Gf));
Tn = Ln/(1+Ln);
Sn = 1/(1+Ln);

L_II = minreal(zpk(inv(S_II)-1));
%nichols(Ln)
omega2 = logspace(-10,10,1000);
[magLn,phaseLn] = nichols(Ln,omega2);
magLn = squeeze(magLn);
phaseLn = squeeze(phaseLn);
plot(phaseLn,20*log10(magLn),'LineWidth',2)

figure("Name","Big Picture");
%superimpose Gc
myngridst(Tp0,Sp0); 
hold on,grid on
[magLn,phaseLn]=nichols(Ln,logspace(-150, 150,3000));
magLn=squeeze(magLn);
phaseLn=squeeze(phaseLn);
plot(phaseLn,20*log10(magLn),'LineWidth',2,'Color','blue')
Lnmod=minreal(Gcmod*Ga*Gpn*Gs*Gf);
[magLnmod,phaseLnmod]=nichols(Lnmod,logspace(-150, 150,3000));
magLnmod=squeeze(magLnmod);
phaseLnmod=squeeze(phaseLnmod);
plot(phaseLnmod,20*log10(magLnmod),'LineWidth',2,'Color','red')
axis([-360 0 -150 150])
title('Problem 1 -- Loop functions (big picture)','Interpreter','latex','FontSize',12);
xlabel('Phase (deg)','Interpreter','latex','FontSize',12);
ylabel('Magnitude (dB)','Interpreter','latex','FontSize',12); 
legend('','','','','$L(G_c)$','$L(G_{c,mod})$','Interpreter','latex','FontSize',12)

%bounds on the loop function L(s)
figure("Name","Bounds on L(s)")
[magLn,~]=bode(Ln,omega);
magLn=squeeze(magLn); 
semilogx(omega,20*log10(magLn),'LineWidth',1);
hold on, grid on

%lower bound
wlb=logspace(-4,-1,100);
[magWs,~]=bode(1/WS_inv,wlb); 
magWs=squeeze(magWs);
[magWu,~]=bode(Wu,wlb); 
magWu=squeeze(magWu);
lb=squeeze((magWs)./(1-magWu));
semilogx(wlb,20*log10(lb),'Color','red','LineWidth',2);

%upper bound
wub=logspace(1,4,100);
[magWs,~]=bode(1/WS_inv,wub); 
magWs=squeeze(magWs);
[magWu,~]=bode(Wu,wub); 
magWu=squeeze(magWu);
ub=squeeze(abs((1-magWs)./(magWu)));
semilogx(wub,20*log10(ub),'Color','green','LineWidth',2);

% check of robust stability ||WuTn||_inf < 1 -> ||Tn||_inf < ||Wu_inv||_inf
[magWu_inv,~]=bode(1/Wu,omega); 
magWu_inv=squeeze(magWu_inv);
figure("Name","robust stability check");
semilogx(omega,20*log10(magWu_inv),'LineWidth',1, 'Color','red');
hold on, grid on
[magTn,~]=bode(Tn,omega);
magTn=squeeze(magTn);
semilogx(omega,20*log10(magTn),'LineWidth',1)
title(['Robust stability ' ...
    '$\Vert W_u T_n \Vert_\infty < 1$'],'Interpreter','latex','FontSize',14);
xlabel('Frequency (rad/s)','Interpreter','latex','FontSize',12);
ylabel('Magnitude (dB)','Interpreter','latex','FontSize',12);

[magWu,~]=bode(Wu,omega); 
magWu=squeeze(magWu);

Hinf_WuTn=norm(abs(magWu.*magTn),inf)

%check of nominal performance ||Sn||_inf < ||WS_inv||_inf and ||Tn||_inf <
%||WT_inv||_inf
figure("Name","Nominal performance check(1)");
[magWT_inv,~]=freqresp(WT_inv,omega); 
magWT_inv=abs(squeeze(magWT_inv)); 
semilogx(omega,20*log10(magWT_inv),'LineWidth',1.5);
hold on, grid on
[magTn,~]=bode(Tn,omega);
magTn=squeeze(magTn);
semilogx(omega,20*log10(magTn),'LineWidth',1)
title(['Nominal performance ' ...
    '$\Vert W_T T_n \Vert_\infty < 1$'],'Interpreter','latex','FontSize',14);
xlabel('Frequency (rad/s)','Interpreter','latex','FontSize',12);
ylabel('Magnitude (dB)','Interpreter','latex','FontSize',12);


[magWT,~]=freqresp(1/WT_inv,omega); 
magWT=abs(squeeze(magWT)); 
Hinf_WtTn=norm(abs(magWT.*magTn),inf)

figure("Name","Nominal performance check(2)");
[magWs,~]=freqresp(WS_inv,omega);
magWs=abs(squeeze(magWs));
semilogx(omega,20*log10(magWs),'LineWidth',1.5);
hold on, grid on
[magSn,~]=bode(Sn,omega);
magSn=squeeze(magSn);
semilogx(omega,20*log10(magSn),'LineWidth',1)
title(['Nominal performance ' ...
    '$\Vert W_S S_n \Vert_\infty < 1$'],'Interpreter','latex','FontSize',14);
xlabel('Frequency (rad/s)','Interpreter','latex','FontSize',12);
ylabel('Magnitude (dB)','Interpreter','latex','FontSize',12);

[magWs,~]=freqresp(1/WS_inv,omega);
magWs=abs(squeeze(magWs));

Hinf_WsSn=norm(abs(magWs.*magSn),inf)

%check of robust performance || |WS*Sn|+|Wu*Tn| ||_inf < 1
[magWuTn,~]=bode(minreal(Wu*Tn),omega); 
magWuTn=squeeze(magWuTn);
Ws=1/WS_inv;
[magWsSn,~]=bode(minreal(Ws*Sn),omega); magWsSn=squeeze(magWsSn);
figure("Name","Check of robust performance"); 
semilogx(omega,20*log10(magWuTn),'LineWidth',1,'LineStyle','--'); 
hold on, grid on
semilogx(omega,20*log10(magWsSn),'LineWidth',1,'LineStyle','--');
sum=magWuTn+magWsSn; 
semilogx(omega,20*log10(sum),'LineWidth',2,'Color','red')
title(['Robust performance ' ...
    '$\Vert \vert W_u T_n \vert + \vert W_S S_n\vert \Vert_\infty< 1$'],'Interpreter','latex','FontSize',14);
xlabel('Frequency (rad/s)','Interpreter','latex','FontSize',12);
ylabel('Magnitude (dB)','Interpreter','latex','FontSize',12);

Hinf_WuTn_WuSn = norm(sum,inf)

%check of transient requirements
tr0 = 3;
ts0 = 12;
figure("Name","Transient Analysis")
step(Tn)
hold on, grid on
xline(tr0), xline(ts0), yline(Kd+Kd*0.1), yline(Kd+0.05*Kd, 'r--'), yline(Kd-0.05*Kd, 'r--'),
yline(0.1*Kd, 'm--'), yline(0.9*Kd, 'm--')
title('Transient analysis','Interpreter','latex','FontSize',14);
xlabel('Time (s)','Interpreter','latex','FontSize',12);
ylabel('$y(t)$','Interpreter','latex','FontSize',12);

%check of (S2)
R0_t=R0; 
Da0_t=0; 
ap_t=0; 
as_t=0; 
Tsim=100*ts0;
out=sim('general_model.slx');
figure("Name","e_r")
plot(out.e_r.time, out.e_r.data);
grid on

%check of (S3)
R0_t=0; 
Da0_t=Da0;
ap_t=0; 
as_t=0; 
Tsim=100*ts0;
out=sim('general_model.slx');
figure("Name","e_da")
plot(out.e_yd.time, out.e_yd.data);
grid on

%check of (S4)
R0_t=0; 
Da0_t=0;
ap_t=ap; 
as_t=0; 
Tsim=200*ts0;
out=sim('general_model.slx');
figure("Name","e_dp")
plot(out.e_yd.time, out.e_yd.data);
grid on

%check of (S5)
R0_t=0; 
Da0_t=0;
ap_t=0; 
as_t=as;
Tsim=5;
out=sim('general_model.slx');
figure("Name","e_ds")
plot(out.e_yd.time, out.e_yd.data);

grid on

% % ----------- mu-analysis ---------------------

Wk = (Ku(2)-Kn)/Kn;
Wp1 = (p1u(2)-p1n)/p1n;
Wp2 = (p2u(2)-p2n)/p2n;

%RP on S3 (polynomial disturbance) - structured
WS_mu = minreal(1/(s*S_star0_S3));
[An,Bn,Cn,Dn] = linmod('N_WS');
Ns = ss(An,Bn,Cn,Dn);
Ns = minreal(zpk(Ns), 1e-3);
[An,Bn,Cn,Dn] = ssdata(Ns); 
N = pck(An,Bn,Cn,Dn);
omega = logspace(-8,-4,500);
Nf = frsp(N, omega);
deltaset = [-1 0; -1 0; -1 0; 1 1];
mubnds = mu(Nf, deltaset);
vplot('liv,m', mubnds);
title('Robust performance for $\omega\to0$', ...
    'Interpreter','latex','FontSize',14);
subtitle('Polynomial disturbance $d_a(t)$ - Structured Uncertainty', ...
    'Interpreter','latex');
legend('upper bound','lower bound','Interpreter','latex');
grid on

%RP on S4 (sinusoidal disturbance) - structured
WS_mu = tf(1,MS_LF);
[An,Bn,Cn,Dn] = linmod('N_WS');
Ns=minreal(zpk(ss(An,Bn,Cn,Dn)));
[An,Bn,Cn,Dn] = ssdata(Ns); 
N = pck(An,Bn,Cn,Dn);
omega = logspace(log10(0.01*wp_plus),log10(wp_plus),500);
Nf = frsp(N,omega);
deltaset = [-1,0; -1,0; -1,0; 1,1];
mubnds = mu(Nf,deltaset);
figure("Name","mu-analysis: RP on S4");
vplot('liv,m', mubnds);
title('Robust performance in the low frequency range', ...
    'Interpreter','latex','FontSize',14);
subtitle('Sinusoidal disturbance $d_p(t)$ - Structured', ...
    'Interpreter','latex');
legend('upper bound','lower bound','Interpreter','latex');
grid on

%RP on S5 (sinusoidal disturbance) - structured
WT_mu = tf(1, MT_HF);
[An,Bn,Cn,Dn] = linmod('N_WT');
Ns = minreal(zpk(ss(An,Bn,Cn,Dn)));
[An,Bn,Cn,Dn] = ssdata(Ns); 
N = pck(An,Bn,Cn,Dn);
omega = logspace(log10(ws_minus),log10(100*ws_minus),500);
Nf = frsp(N,omega);
deltaset = [-1,0; -1,0; -1,0; 1,1];
mubnds = mu(Nf,deltaset);
figure("Name","mu-analysis: RP on S5");
vplot('liv,m', mubnds);
title('Robust performance in the high frequency range', ...
    'Interpreter','latex','FontSize',14);
subtitle('Sinusoidal disturbance $d_s(t)$ - Structured','Interpreter','latex');
legend('upper bound','lower bound','Interpreter','latex');
grid on

%Stability check - Structured
[An,Bn,Cn,Dn] = linmod('N11');
Ns = minreal(zpk(ss(An,Bn,Cn,Dn)));
[An,Bn,Cn,Dn]=ssdata(Ns); 
N = pck(An,Bn,Cn,Dn);
omega = logspace(-3,3,500);
Nf = frsp(N,omega);
deltaset = [-1,0; -1,0; -1,0];
mubnds = mu(Nf,deltaset);
figure("Name","mu-analysis: RS");
vplot('liv,m', mubnds);
title('Robust stability test', ...
    'Interpreter','latex','FontSize',14);
subtitle('Structured uncertainty', ...
    'Interpreter','latex');
legend('upper bound','lower bound','Interpreter','latex');
grid on












































