% EES-40 2023-2 Lab1 semana 4 Jacques Waldmann 18 agosto 2023
% realização diagonal por blocos: 10/(s+10) é dinâmica de primeira ordem em cascata com
% transformação para base modal da dinâmica de segunda ordem
% G(s)= 10/(s+10) * 100/(s+0,5+10j)/(s+0,5-10j)
% o grupo deve usar outra realização
clear all;clc;close all
A=[-10 0 0;-10 -.5 10;0 -10 -.5];
B=[10;0;0];
C=[0 0 1];

[Gtfnum Gtfden]=ss2tf(A,B,C,0);
Gtf=tf(Gtfnum,Gtfden);
step(Gtf)

[An,Bn,Cn,Dn] = tf2ss(Gtfnum,Gtfden)

[Gtfnumn Gtfdenn]=ss2tf(An,Bn,Cn,0);
Gtfn=tf(Gtfnumn,Gtfdenn);
step(Gtfn)
print('figs/Gplanta.eps', '-depsc');
%% augmented dynamics with n=4 to model integral control action
% augmenting state vector x with component x_i: time integral of r-y (output 
% error with respect to the reference) - augmented state vector is xaug augmented 
% dynamics: xaug_dot=Aaug*xaug+Baug*u+Gaug*r y=Caug*xaug

Aaug=[An zeros(3,1);-Cn 0]
Baug=[Bn;0] % % control u input matrix
Gaug=[zeros(3,1); 1] % reference r input matrix
Caug=[Cn 0]
%% ITAE criterion n=4 poles, w0 is the desired bandwidth [rd/s]
% desired poles are (s/w0+.424+-1.263j) and (s/w0+.626+-.414j) Ogata segunda 
% edição pg. 240 Tabela 4.2 vide também arquivo PDF da KNTU a esse respeito

w0=1*3*pi; % rd/s (Hz->rd/s) selected here is 1 Hz
X=['desired bandwidth ',num2str(w0),' rd/s'];
disp(X);
s=tf('s');
polinmf=[1 2.1*w0 3.4*w0^2 2.7*w0^3 w0^4];
p=roots(polinmf); % 4 polos projetados para malha fechada via realimentação de estado aumentado
X=['vetor de estado aumentado: polos de malha fechada desejados ',num2str(p(1)),'  ',num2str(p(2)),'  ',num2str(p(3)),'  ',num2str(p(4))];
disp(X);
%%
Kaug=place(Aaug,Baug,p);
K=Kaug(1:3); % feedback gain to be used with negative feedback => u= - Kaug*xaug includes integral action
Ki=-Kaug(4); % gain Ki weighs the integral of tracking error r.1(t)-y
X=['augmented state feedback gain vector Kaug ',num2str(Kaug)];
disp(X);
Aaugmf=Aaug-Baug*Kaug;
pmf=eig(Aaugmf);
X=['polos de malha fechada obtidos ',num2str(pmf(1)),'  ',num2str(pmf(2)),'  ',num2str(pmf(3)),'  ',num2str(pmf(4))];
disp(X);
%% closed-loop response in steady state

step_mag=5; % magnitude que satura amp-op
M=[An Bn;-Cn 0];
xureg=inv(M)*[zeros(3,1);-step_mag] % plant model steady state is xureg(1:3) and corresponding control is xureg(4)
xireg=1/Ki*(xureg(4)+K*xureg(1:3)) % augmented state component xi in steady state
%% closed-loop response

SysMFpos=ss(Aaugmf,Gaug,Caug,0); % output is the third component of the augmented state vector
StepInfo=stepinfo(SysMFpos,'RiseTimeThreshold',[0 1]);
X=['Overshoot [%] ',num2str(StepInfo.Overshoot),' Rise time 0-100% [s] ',num2str(StepInfo.RiseTime)];
disp(X);

opt = stepDataOptions('StepAmplitude',step_mag);
[ystep,tstep,xstep]=step(SysMFpos,opt);

figure; % first state vector xaug component
plot(tstep,xstep(:,1));
title(['xaug(1) em malha fechada - r(t)= ',num2str(step_mag),'*1(t)']);
xlabel('segundos');
grid;
print('figs/x1mf.eps', '-depsc');

figure; % second state vector xaug component
plot(tstep,xstep(:,2));
title(['xaug(2) em malha fechada - r(t)= ',num2str(step_mag),'*1(t)']);
xlabel('segundos');
grid;
print('figs/x2mf.eps', '-depsc');

figure; % third state vector xaug component
plot(tstep,xstep(:,3));
title(['xaug(3) em malha fechada - r(t)= ',num2str(step_mag),'*1(t)']);
xlabel('segundos');
grid;
print('figs/x3mf.eps', '-depsc');

figure; % fourth state vector xaug component is integral of r-y
plot(tstep,xstep(:,4));
title(['xaug(4) em malha fechada - r(t)= ',num2str(step_mag),' *1(t)']);
xlabel('segundos');
grid;
print('figs/x4mf.eps', '-depsc');

figure; % sinal de saída em malha fechada
plot(tstep,ystep);
title(['saída em malha fechada - vide matriz Caug - r(t)= ',num2str(step_mag),' *1(t)']);
xlabel('segundos');
grid;
print('figs/saida.eps', '-depsc');

figure;   % sinal de controle em malha fechada com referência degrau
% control law: u=-Kaug(1:3)*xstep(1:3,:)'+Ki*x_step(4,:)'=-Kaug*xstep'
% notice that xstep(4,:)=integral of(r-y) in closed-loop simulation
xstep=xstep'; % estado em malha fechada
cntrl=-Kaug*xstep;
plot(tstep,cntrl);
title(['controle[V] em malha fechada - r(t)= ',num2str(step_mag),' *1(t)']);
xlabel('segundos');
grid;
print('figs/sinaldecontrole.eps', '-depsc');
%% stability margins - dictated by the open loop L
% Notice that u=Kaug*xaug is output and negative feedback is assumed

[Num,Den]=ss2tf(Aaug,Baug,Kaug,0); % nominal plant model augmented with integral action
L=tf(Num,Den);                     % is the open loop tf L with u=Kaug*xaug as output and negative feedback is
% assumed in the stability margins analysis of the open loop L
figure;
nyquistplot(L); % use zoom in critical point in the Nyquist plot and show margins
title('Nyquist plot - open loop L');
figure;
margin(L);
% Design is not robust because nominal plant model L comes too close to -1,0j and the closed loop
% can become unstable due to uncertainty in L
D=1+L; % distance from critical point in complex plane (-1+0j) to open loop tf L
S=1/D; % sensitivity tf
figure;
bodemag(S);
title('Magnitude da função de sensibilidade |S(j\omega)|[dB]'); % projeto não é robusto
figure;
bodemag(D);
title('função distância ao ponto crítico |D(j\omega|[dB]');
