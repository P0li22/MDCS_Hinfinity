clear all
close all
clc

%% --------------nominal plant----------------
s = tf('s');
Ku = [80; 120]; p1u = [0.8; 1.2]; p2u = [2; 7];
Kn = (Ku(1)+Ku(2))/2; p1n = (p1u(1)+p1u(2))/2; p2n = (p2u(1)+p2u(2))/2;
Gpn = minreal(zpk(Kn/(4.5*s*(1+s/p1n)*(1+s/p2n))));
Kp = dcgain(s*Gpn);

%% ---------------requirements--------------------
%S1
Ga = 0.112; Gs = 0.5; 
Kd = 8;
Gf = 1/(Gs*Kd);

%S2
R0 = 1;
nu_s2 = 0;
rho_r = 1.5e-1;
Sstar0_S2 = rho_r / abs(R0*Kd);

%S3
rho_a = 2.14;
Da0 = 1.5e-3;
nu_s3 = 1;
Sstar0_s3 = rho_a/abs(Kp*Da0);

%S4
rho_p = 5.1e-3;
ap = 16e-2;
wp = 0.03;
MS_LF = rho_p/ap;
MS_LF_dB = 20*log10(MS_LF);
wL = wp/sqrt(MS_LF);

%S5
rho_s = 1.6e-3;
as = 2e-1;
ws = 60;
MT_HF = rho_s*Gs/as;
MT_HF_dB = 20*log10(MT_HF);
wH = ws*sqrt(MT_HF);

%S6-S8
s_hat = 0.13;
zetaII = abs(log(s_hat))/sqrt(pi^2+(log(s_hat))^2);
Sp0 = 2*zetaII*sqrt(2+4*zetaII^2+2*sqrt(1+8*zetaII^2))/(sqrt(1+8*zetaII^2)+4*zetaII^2-1);
Sp0_dB = 20*log10(Sp0);
Tp0 = 1/(2*zetaII*sqrt(1-zetaII^2));
Tp0_dB = 20*log10(Tp0);

tr0 = 1.8;
wnII1 = (pi-acos(zetaII))/(tr0*sqrt(1-zetaII^2));

ts0 = 6; alpha = 0.05;
wnII2 = (-log(alpha))/(ts0*zetaII);

wnII = max(wnII1, wnII2);
wc = wnII*sqrt(sqrt(1+4*zetaII^4)-2*zetaII^2);

%% ----------------------Ws----------------------
nu = 1; p = 1;
omega = logspace(-3, 3, 1000);
Sstar0 = Sstar0_s3;
Sto0 = minreal(s^(nu+p)*Sstar0);
[magSto0, ~] = bode(Sto0, omega); magSto0 = squeeze(magSto0);
SII = minreal(s*(s+2*zetaII*wnII)/(s^2+2*zetaII*wnII*s+wnII^2));
[magSII, ~] = bode(SII, omega); magSII = squeeze(magSII);

a = Sstar0;
p1 = 0.012;
zetaws = 0.78;
w1 = 0.79*wnII;
z1 = (a*p1*w1^2)/Sp0;
WsInv = minreal( (a*s^(nu+p)*(1+s/z1))/((1+s/p1)*(1+s*(2*zetaws/w1)+s^2/w1^2)) );
[magWsInv, ~] = bode(WsInv, omega); magWsInv = squeeze(magWsInv);

figure("Name", "Ws")
semilogx(omega, 20*log10(magSto0), 'r--');
hold on, grid on
semilogx(omega, 20*log10(magSII), 'b--');
semilogx(omega, 20*log10(magWsInv), 'LineWidth',2,'Color','k');
xline(wp), yline(MS_LF_dB), yline(Sp0_dB, 'm')

Ws = minreal(zpk(inv(WsInv)))
[magWs, ~] = bode(Ws, omega); magWs = squeeze(magWs);

%% -----------------Wt---------------------
TII = minreal(1 / (1+s*(2*zetaII/wnII)+s^2/wnII^2));
[magTII, ~] = bode(TII, omega); magTII = squeeze(magTII);

wT = ws*10^((MT_HF_dB-Tp0_dB)/40)
WtInv = minreal(Tp0 / (1+s*(1.414/wT)+s^2/wT^2));
[magWtInv, ~] = bode(WtInv, omega); magWtInv = squeeze(magWtInv);

figure("Name", "Wt")
semilogx(omega, 20*log10(magTII), 'r--');
hold on, grid on
semilogx(omega, 20*log10(magWtInv), 'LineWidth',2,'Color','k');
xline(ws), yline(MT_HF_dB), yline(Tp0_dB, 'm')

Wt = minreal(zpk(inv(WtInv)))
[magWt, ~] = bode(Wt, omega); magWt = squeeze(magWt);

%% ---------------------Wu---------------------------
% nsamp = 10;
% gridK = linspace(Ku(1), Ku(2), nsamp);
% gridp1 = linspace(p1u(1), p1u(2), nsamp);
% gridp2 = linspace(p2u(1), p2u(2), nsamp);
% matDelta = zeros(nsamp^3, 1000); cnt = 1;
% 
% %cloud of TFs
% figure()
% for ii = 1:nsamp
%     for jj = 1:nsamp
%         for kk = 1:nsamp
%             K_ = gridK(ii);
%             p1_ = gridp1(jj);
%             p2_ = gridp2(kk);
%             Gp_ = minreal(zpk(K_/(4.5*s*(1+s/p1_)*(1+s/p2_))));
%             Delta = (Gp_/Gpn)-1;
%             [magDelta, ~] = bode(Delta, omega); 
%             magDelta = squeeze(magDelta);
%             matDelta(cnt, :) = magDelta';
%             cnt = cnt+1;
%             semilogx(omega, 20*log10(magDelta));
%             hold on
%         end
%     end
% end
% 
% maxDelta = max(matDelta)';
% WuMag = vpck(maxDelta, omega);
% fit = fitmag(WuMag);
% [A, B, C, D] = unpck(fit);
% Wu = ss(A, B, C, D);
% Wu = minreal(zpk(Wu))
% 
% figure("Name", "Wu")
% for ii = 1:nsamp
%     for jj = 1:nsamp
%         for kk = 1:nsamp
%             K_ = gridK(ii);
%             p1_ = gridp1(jj);
%             p2_ = gridp2(kk);
%             Gp_ = minreal(zpk(K_/(4.5*s*(1+s/p1_)*(1+s/p2_))));
%             Delta = (Gp_/Gpn)-1;
%             [magDelta, ~] = bode(Delta, omega); 
%             magDelta = squeeze(magDelta);
%             matDelta(cnt, :) = magDelta';
%             cnt = cnt+1;
%             semilogx(omega, 20*log10(magDelta));
%             hold on
%         end
%     end
% end
% semilogx(omega, 20*log10(magWu), 'LineWidth',2,'Color','k');

Wu = minreal(zpk( 1.2391*(s+0.3358)*(s^2 + 5.999*s + 11.92)/((s+5.867)*(s+4.003)*(s+1.053))))
[magWu, ~] = bode(Wu, omega); magWu = squeeze(magWu);

%% --------------------------Hinf LMI-----------------------------

%W1mod
W1 = Ws;
wc = 1;
lambda = 0.01*wc;
W1mod = minreal(W1*(s/(s+lambda))^(nu+p));

%W2mod
figure("Name","Wu vs Wt");
bodemag(Wt), hold on, bodemag(Wu), grid on
W2 = Wt;
W2mod = tf(1, Tp0);

% gen plant RS NP
[A, B, C, D] = linmod("P7_genPlant_RS_NP");
M = ltisys(A, B, C, D);

%insert zeros
M = sderiv(M, 2, [1/wT 1]);
M = sderiv(M, 2, [1/wT 1]);

%LMI
[gopt, Gcmod] = hinflmi(M, [1 1], 0, 0.01, [0 0 0]);
[A, B, C, D] = ltiss(Gcmod);
Gcmod = ss(A, B, C, D);
Gcmod = minreal(zpk(Gcmod), 1e-4)
Gc = minreal(zpk( Gcmod*((s+0.01005)/s)*((s+0.008179)/(s+4.045e-06))*(1+s/1.15e05) ), 1e-3)

%% ------------------check of performance------------------------
% Nichols plot -> NS
omega = logspace(-100, 100, 1000);
Lnmod = minreal(Gcmod*Gpn*Ga*Gs*Gf);
Ln = minreal(Gc*Gpn*Ga*Gs*Gf);
[magLnmod, pLnmod] = bode(Lnmod, omega); 
magLnmod = squeeze(magLnmod); pLnmod = squeeze(pLnmod);
[magLn, pLn] = bode(Ln, omega); magLn = squeeze(magLn); pLn = squeeze(pLn);

figure("Name", "Nichols")
myngridst(Tp0, Sp0);
hold on, grid on
semilogx(pLnmod, 20*log10(magLnmod), 'LineWidth',2,'Color','r');
semilogx(pLn, 20*log10(magLn), 'LineWidth',2,'Color','k');

Tn = minreal(Ln/(1+Ln));
Sn = minreal(1/(1+Ln));

% RS 
omega = logspace(-3, 3, 1000);
figure("Name", "RS")
bodemag(inv(Wu), 'r');
hold on
bodemag(Tn, 'k');
grid on

Hinf_WuTn = norm(minreal(Wu*Tn), inf)

% NP 
figure("Name", "NP Sn");
bodemag(WsInv, 'r');
hold on
bodemag(Sn, 'k');

Hinf_WsSn = norm(minreal(Ws*Sn), inf)

figure("Name", "NP Tn");
bodemag(WtInv, 'r');
hold on
bodemag(Tn, 'k');

Hinf_WtTn = norm(minreal(Wt*Tn), inf)

%time domain
figure("Name", "transient analysis")
step(Tn);
yline(1.13), yline(0.95, 'b'), yline(1.05, 'b'), yline(Kd, 'r')
xline(tr0, 'm'), xline(ts0, 'm--')

%% ---------------mu analysis----------------------------
Wk = (Ku(2)-Kn)/Kn;
Wp1 = (p1u(2)-p1n)/p1n;
Wp2 = (p2u(2)-p2n)/p2n;

%RP S3
Ws_mu = minreal(1/(s^(nu_s3+p)*Sstar0_s3));
[A, B, C, D] = linmod("N_WS");
Ns = ss(A, B, C, D);
Ns = minreal(Ns);
[A, B, C, D] = ssdata(Ns);
N = pck(A, B, C, D);
omega = logspace(-8, -3, 500);
Nf = frsp(N, omega);
deltaset = [-1 0; -1 0; -1 0; 1 1];
mubounds = mu(Nf, deltaset);
figure()
vplot('liv,m', mubounds);
title("RP check S3");
grid on

%RP S4
Ws_mu = tf(1, MS_LF);
[A, B, C, D] = linmod("N_WS");
Ns = ss(A, B, C, D);
Ns = minreal(Ns);
[A, B, C, D] = ssdata(Ns);
N = pck(A, B, C, D);
omega = logspace(log10(0.01*wp), log10(wp), 500);
Nf = frsp(N, omega);
deltaset = [-1 0; -1 0; -1 0; 1 1];
mubounds = mu(Nf, deltaset);
figure()
vplot('liv,m', mubounds);
title("RP check S4");
grid on

%RS
[A, B, C, D] = linmod("N11");
Ns = ss(A, B, C, D);
Ns = minreal(Ns);
[A, B, C, D] = ssdata(Ns);
N = pck(A, B, C, D);
omega = logspace(-3, 3, 500);
Nf = frsp(N, omega);
deltaset = [-1 0; -1 0; -1 0];
mubounds = mu(Nf, deltaset);
figure()
vplot('liv,m', mubounds);
title("RS check");
grid on