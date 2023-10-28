% EES-40 2023-2 Lab3 semana 12 Jacques Waldmann outubro 2023
% Controle ótimo LQG de sistema SISO 
% Polos de malha fechada selecionados via root locus simétrico (LGR-S) 
% Robustez do controle ótimo LQR - função de sensibilidade e as margens de estabilidade
% Solução LQG - combina LQR com realimentacao do estado observado pelo 
% estimador LQE (filtro de Kalman em regime) 
% Função de sensibilidade do controle ótimo LQG - degradação das margens de estabilidade
% Desempenho LQG e o compromisso com o esforço de controle e a passagem de ruído de medida 

% Realização em cascata de primeira ordem com  modal de segunda ordem
% G(s)= 10/(s+10) * 100/(s+0,5+10j)/(s+0,5-10j) = Y(s)/U(s) é a malha aberta
clear all;clc;close all
% cada grupo com realização {A,B,C} distinta
A=[-10 0 0;-10 -.5 10;0 -10 -.5];
B=[10;0;0];
C=[0 0 1];
Gss=ss(A,B,C,0);
[Gtfnum Gtfden]=ss2tf(A,B,C,0);
Gtf=tf(Gtfnum,Gtfden); % malha aberta, sem ação integral
figure;
step(Gtf)

T = [[1, 2, 0]; [0, 2, 0]; [0, 0, 1]]
det(T)
A = inv(T) * A * T
B = inv(T) * B
C = C * T
Gss=ss(A,B,C,0);
[Gtfnum Gtfden]=ss2tf(A,B,C,0);
Gtf=tf(Gtfnum,Gtfden); % malha aberta, sem ação integral
figure;
step(Gtf)

%% augmented state xaug dynamics 
% includes xi: time integral of r-y (output error wrt to reference) 
% xaug_dot=Aaug*xaug+Baug*u+Gaug*r
% y=Caug*xaug
Aaug=[A zeros(3,1);-C 0];
Baug=[B;0];
Gaug=[zeros(3,1); 1]; % reference r input matrix  
Caug=[C 0]; % physical output y
% C_xi is the infinite-horizon time integral of the output tracking error 
% for use in the quadratic cost function J_LQR
C_xi=[0 0 0 1];  
s=tf('s'); % Laplace transform - variable s for later use in transfer functions  
%% lugar geométrico das raízes simétrico (LGR-S)
% symmetric root locus (SRL) LGR-S to select closed-loop poles according to rho^(-1)
% minimize J_LQR=integral (xi^2 + rho u^2)dt from t=0 to +inf
[numAug,denAug]=ss2tf(Aaug,Baug,C_xi,0);
Gaugtf=tf(numAug,denAug); % malha aberta, inclui ação integral, com "saída" xi ponderada em J_LQR
% lqr function yields Klqr(1:4)=[K  -Ki] such that (be aware of the negative sign)
%      u= - Klqr * [x'  xi']' is the control signal optimizing J_LQR 
% while driving the augmented linear dynamics and tracking the reference.  
% LQR design assumes the actual state is available 
% Closed-loop pole selection according to symmetric root locus with output xi in J_LQR  
num1=numAug; % numerador de Gaugtf(-s)
den1=denAug.*[1 -1 1 -1 1]; % denominador de Gaugtf(-s), que tem xi como "saída"
num_s=conv(numAug,num1); % numerador de Gaugtf(s)*Gaugtf(-s)
den_s=conv(denAug,den1); % denominador de Gaugtf(s)*Gaugtf(-s)
sys_s=tf(num_s,den_s); % open loop tf Gaugtf(s)*Gaugtf(-s) for SRL
rlocus(sys_s); % choose with mouse acceptable gain from the stable closed-loop poles in SRL asymptotes
%%
%GainLGR_S=input(' selected SRL gain (inverse of rho seen in J-cost functional ');
GainLGR_S = 30
rho=1/GainLGR_S
Klqr=lqr(Aaug,Baug,C_xi'*C_xi,rho); % lqr com peso em xi e rho=-1/GainLGR_S no controle, sendo GainLGR_S obtido com LGR simétrico para atender requisitos
%% closing the loop with LQR optimal state feedback gain
Ki=-Klqr(4);
K=Klqr(1:3);
Aaugmf=[A-B*K B*Ki;-C 0]; % closed-loop dynamics matrix
pmf=eig(Aaugmf); % polos de malha fechada
X=['state feedback gain vector K ',num2str(K(1:3)),' output tracking error integral gain Ki ',num2str(Ki)];
disp(X);
X=['attained closed-loop poles ',num2str(pmf(1)),'  ',num2str(pmf(2)),'  ',num2str(pmf(3)),'  ',num2str(pmf(4))];
disp(X);
%% model-based steady state analysis: closed-loop response in steady state 
step_mag=5; % reference step magnitude [V] 
X=['ref magnitude ',num2str(step_mag),' [V]'];disp(X);
M=[A B;-C 0];
xureg=inv(M)*[zeros(3,1);-step_mag] % plant model steady state is xureg(1:3) and control xureg(4)
xireg=1/Ki*(xureg(4)+K*xureg(1:3)) % augmented state component xi in steady state
%% closed-loop, augmented-state feedback simulation of the LQR design
% notice: input is reference r, output are y and xi=xaug(4) 
SysMFpos=ss(Aaugmf,Gaug,[Caug;C_xi],0); 
damp(SysMFpos)  % closed-loop poles natural frequency and damping ratio
opt = stepDataOptions('StepAmplitude',step_mag);
tfinal=2; % final time tf=2[s]
% xstep is the timeline of xaug; ystep is the timeline of xi=xaug(4)
[ystep,tstep,xstep]=step(SysMFpos,tfinal,opt);  

figure; % sinal xaug1(t) em malha fechada
plot(tstep,xstep(:,1));
title(['xaug(1) em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

figure; % sinal xaug2(t) em malha fechada
plot(tstep,xstep(:,2));
title(['xaug(2) em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

figure; % sinal xaug(3) em malha fechada
plot(tstep,xstep(:,3));
title(['xaug(3) em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

figure; % output signal: xaug(4)=xi is the closed-loop time-integral of the output tracking error 
plot(tstep,ystep(:,1:2));
title(['saída y e "saída" xi em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V.s]']);
legend('y','x_I');
xlabel('segundos');
grid;

figure;   % sinal de controle em malha fechada 
xstep_transp=xstep'; % linha temporal do estado aumentado em malha fechada é aqui transposto
cntrl=-Klqr*xstep_transp;
plot(tstep,cntrl);
title(['controle LQR [V] em malha fechada - degrau_{ref}= ',num2str(step_mag),' [V]']);
xlabel('segundos');
grid;

%% LQR design stability margins with augmented state xaug feedback is dictated by 
% the open-loop transfer function tfL(jw) where 
% u=Klqr*xaug is the output and negative feedback is assumed
% projeto LQR,além de estabilidade exponencial assegurada, 
% também garante que tfL(jw) passa fora do círculo unitátio centrado no ponto crítico -1+0*j
% margem de ganho inferior menor que -6 dB, superior de +infinito, e margem de fase superior a +60 graus ou inferior a -60 graus 
% inspecionar quão próximo do ponto crítico passa o modelo nominal 
SysMA=ss(Aaug,Baug,Klqr,0);
figure;
nyquistplot(SysMA); % use zoom in critical point in the Nyquist plot and see gain. phase margins
title('Nyquist plot - open loop tfL(j\omega)');
figure;
margin(SysMA);
% para inspecionar robustez ver se tfL(jw) se aproxima muito do ponto crítico - pode 
% se tornar instável com modelo incerto e o sintoma é envolver o ponto crítico 
[Num,Den]=ss2tf(Aaug,Baug,Klqr,0); 
tfL=tf(Num,Den); % open loop tf with u=Klqr*xaug as output and negative feedback is assumed 
Dlqr=1+tfL; % distance from critical point in complex plane (-1+0j) to open loop tfL(jw)
Slqr=1/Dlqr; % sensitivity tf
figure;
bodemag(Slqr);
title('Magnitude da função de sensibilidade |S_{lqr}(j\omega)|[dB]'); 
figure;
bodemag(Dlqr);
title('função distância ao ponto crítico |D_{lqr}(j\omega|[dB]'); 

%% Lab3 - The LQE observer estimates the plant nominal model state vector 
% Consequently, the observer estimates just the three initial components (x) of the augmented state vector xaug
% The available noisy measurement is 
% y = C*x + v = [0 0 1]*x + v with the suggested realization
% and scalar Gaussian white noise v with PSD R. 

% Augmented nominal plant model with Gaussian white noise uncertainty model is 
% xdot=A*x+B*uhat+F*w
% x_idot=r-y
% w is the scalar white noise in the dynamics model with PSD Q. 
% We use F=B as the model noise w input matrix.
% uhat=-K*xhat+K_i*x_i
% Each group should use A,B,C accordingly consistent with the chosen realization. 

% observer poles are located according to PSDs of Gaussian white, zero-mean model noise w and sensor noise v.
% when Q and R are unknown, common use is to have Q/R as a tuning parameter in SISO models
Q=1; %[PSD units]
% sensor additive scalar Gaussian white noise v with variance R
%R=input('sensor additive white noise PSD ');    % [PSD units]
R = 0.1
[Lo,P,eig_Obs]=lqe(A,B,C,Q,R);

X=['Observer eigenvalues A-Lo*C '];disp(X);
X=num2str(eig_Obs);disp(X);
X=['ganhos do observador Lo'];disp(X);
X=num2str(Lo);disp(X);
%% stability margins - feedback of plant model observed state + integral action on output tracking error r-y 
% observer estimates the augmented plant model state x and xi
% uhat(s)=-(Hy(s)*Y(s)+Hr(s)*R(s)) is the plant input (neglecting input noise to plant)
% negative feedback and feedforward required as shown by the negative sign
% reference does not affect open-loop stability analysis 

% controller
Amf_Obs=[A-Lo*C-B*K  B*Ki;zeros(1,3) 0];
Lo_aug=[Lo;-1];
Hyss=ss(Amf_Obs,Lo_aug,Klqr,0) % feedback path controller; negative feedback required
Hrss=ss(Amf_Obs,Gaug,Klqr,0) % feedforward path controller; negative feedforward required
X=['Continuous-time controller realization - inputs [y;r] and output -uhat with integral action; ***use with negative feedback***'];
disp(X);

% controller state space form concatenates Hyss and Hrss
% output Klqr signal inversion required for use in Simulink and NI-board 
Controlss=ss(Amf_Obs,[Lo_aug Gaug],Klqr,zeros(1,2));

% open-loop transfer function LL_ObsK(jw) stability margins analysis - negative feedback is assumed 
LL_ObsK=minreal(series(Hyss,Gss));
OpenLoop=zpk(LL_ObsK) % poles and zeros of LL_ObsK(jw)
figure;
nyquistplot(LL_ObsK); % use zoom in critical point in the Nyquist plot and show margins
% Notice the 0dB crossover frequency for Tustin prewarping
title('curva de Nyquist - malha aberta LL\_ObsK(j\omega)');
figure;
margin(LL_ObsK); 
legend('LL_ObsK(j\omega) -  see 0dB crossover frequency');

D_ObsK=1+LL_ObsK; % distance from critical point to open loop tf
S_ObsK=1/D_ObsK; % sensitivity tf
figure;
bodemag(S_ObsK);
title('Magnitude da função de sensibilidade |S\_ObsK(j\omega)|[dB]');
figure;
bodemag(D_ObsK);
title('função distância ao ponto crítico |D\_ObsK(j\omega|[dB]')
%% closed-loop Y(s)/R(s) frequency response - bandwidth inspection (-3dB) for 
% sampling period selection  
% inspection of Bode of YOverR reveals closed-loop bandwidth eqs 19 e 20:
A_YOverR=[A-Lo*C zeros(3,3) zeros(3,1);B*K (A-B*K) B*Ki;zeros(1,3) -C 0];
B_YOverR=[zeros(6,1);1]; % r is the input to the closed loop tf
Cymf_YOverR=[zeros(1,5) 1 0]; % output signal y
Cumf_UOverR=-[-K Klqr]; % control signal uhat
YOverRss=ss(A_YOverR,B_YOverR,Cymf_YOverR,0); % closed loop realization [xtilde' x' x_i]' Y/R
UhatOverRss=ss(A_YOverR,B_YOverR,Cumf_UOverR,0); % closed loop realization [xtilde' x' x_i]' U/R

figure;
bode(YOverRss);
title('LQG closed loop bode Y(j\omega)/R(j\omega)- inspect bandwidth frequency at -3 dB crossover');
figure;
step(YOverRss,tfinal,opt); % closed-loop signal y
title('LQG closed loop Y/R - output due to step reference')
figure;
step(UhatOverRss,tfinal,opt); % closed-loop signal y
title('LQG closed loop U/R - control due to step reference')
%% continuous-time closed-loop transfer functions Y/R , U/R - alternative 
% derivation, now with transfer functions
[HytfNum,HytfDen]=tfdata(Hyss); % output of tfdata is in cell representation u/y
Hytf=zpk(tf(HytfNum{1,1},HytfDen{1,1}));
HytfNumMat=cell2mat(HytfNum); % use in continuous-time Simulink file
HytfDenMat=cell2mat(HytfDen);

[HrtfNum,HrtfDen]=tfdata(Hrss); % output of tfdata is in cell representation u/r
Hrtf=zpk(tf(HrtfNum{1,1},HrtfDen{1,1}));
HrtfNumMat=cell2mat(HrtfNum); % use in continuous-time Simulink file
HrtfDenMat=cell2mat(HrtfDen);

YOverR=minreal(series(-Hrtf,feedback(Gtf,Hytf,-1))); % closed loop tf Y/R built with tfs
UOverR=minreal(series(-Hrtf,feedback(1,series(Gtf,Hytf),-1))); % closed loop control tf U/R
%% Verifique resultados no arquivo simulink do sistema continuo no tempo
% condições iniciais distintas entre modelo da planta e observador 
x0=0.1*ones(3,1);
Obs_x0=zeros(4,1);
% Verifique resultados nos scopes do arquivo simulink - sistema continuo no tempo
% antes de prosseguir
%% Tustin discretization with prewarping frequency read at open-loop crossover frequency wc value [rd/s]
%  closed-loop continuous-time bandwidth read in closed-loop Bode -3 dB frequency wb[rd/s] -> transform to [Hz]
wb=input('-3 dB bandwidth frequency [rd/s] from Bode closed loop Y(j\omega)/R(j\omega) ');
fb=wb/2/pi; % Hz
% NI-board maximum D/A sampling frequency is 250Hz
Ts=1/(20*fb); % sufficiently fast sampling time [s] to mitigate Tustin-induced oscillations
save Ts.mat Ts % for use with Tustin discretization and in NI board driver script
%  open-loop continuous-time crossover frequency wc[rd/s] read in open-loop Bode of Y/-Uhat 
wc=input('0 dB crossover frequency [rd/s] from Bode open loop -Uhat(j*\omega)/Y(j*\omega)*G(j*\omega) ');
optcd=c2dOptions('Method', 'tustin', 'PrewarpFrequency',wc);
% Controller C(s) - transfer matrix with feedback and feedforward transfer functions
% C(s) = [Hytf(s) Hrtf(s)]
% Input: [y;r] 2 inputs
% Output: uhat= -Hytf(s)*Y(s)-Hrtf(s)*R(s) ***negative signs required***
% Tustin discretization with prewarping for use in simulink and with the bench plant and NI-board 
Hyd=c2d(Hytf,Ts,optcd) % discrete feedback tf
Hrd=c2d(Hrtf,Ts,optcd) % discrete feedforward tf
[HydNum,HydDen]=tfdata(Hyd); % output of tfdata is in cell representation u/y
HydNumMat=cell2mat(HydNum); % use in continuous-time Simulink file
HydDenMat=cell2mat(HydDen);

[HrdNum,HrdDen]=tfdata(Hrd); % output of tfdata is in cell representation u/r
HrdNumMat=cell2mat(HrdNum); % use in continuous-time Simulink file
HrdDenMat=cell2mat(HrdDen);

% Discretized state space representation of continuous-time controller.
% Controlss embeds the inverted sign - ***ready for use in NI-board drive script***
% reminder: Controlss=ss(Amf_Obs,[Lo_aug Gaug],Klqr,zeros(1,2));
Controlss_d=c2d(-Controlss,Ts,optcd)  % discrete controller state space representation  
save Controlss_d.mat Controlss_d      % for use in NI-board drive script

% continuous-time controller+observer; output: uhat=-Klqr*xaughat and xaughat
% negative feedback is embedded in its observation matrix 
% for use in Simulink file
ObsControlss=ss(Amf_Obs,[Lo_aug Gaug],[-Klqr;eye(4)],zeros(5,2)); 
%ObsControlss=canon(ObsControlss,'modal',CONDT) % modal state space (ss) representation
% discretized controller+observer state space representation 
% *** ready for use in the sampled-data control Simulink file ***
ObsControlss_d=c2d(ObsControlss,Ts,optcd) 
ObsControlss_dA=ObsControlss_d.A;
ObsControlss_dB=ObsControlss_d.B;
ObsControlss_dC=ObsControlss_d.C;
ObsControlss_dD=ObsControlss_d.D;
%%
% 1. rodar arquivo simulink do sistema amostrado e verificar resultados nos
% scopes
% 2. rodar script da placa NI para controle discreto no tempo da planta na
% bancada e verificar resultados
