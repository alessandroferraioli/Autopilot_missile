%% Clear all memory , close all plot and clear the command

clear all;
close all;
clc;
%% Load model and variations


load('model.mat');
load('variations.mat');
load('lqr_ltr.mat');
syms s;
LTR_control=control_2;

la_temp=top_dA;
lm_temp=top_dMout;



%% Define nominal matrixs

An=nominal_model.A
Bn=nominal_model.B
Cn=nominal_model.C
Dn=nominal_model.D

%% check 
disp(ctrb(An,Bn))
disp(rank(ctrb(An,Bn)))



disp(ctrb(An',Cn'))
disp(rank(ctrb(An',Cn')))

%% Check of zero trasmission
tzero(nominal_model)

% valori singolari 
figure
sigma(tf_nominal_model,omega)
grid on
title('Singular values of nominal model');

%% gamma values

gamma1=1;%on S0
gamma2=0.25;% on U0
gamma3=1;%on T0

la=gamma2*la_temp;
lm=gamma3*lm_temp;

%% Calculate S0

close all;
clc;
F0 = series(LTR_control,nominal_model);

%Chiudo il feedback( con il - di default tra F0 ed anello in retroazione
%istantanea--> y=F0*u , ma u=-y ( il - del feedback) . Allora banalmente
%segue che y = (I+F0)^-1 * u , con appunto (I+F0)^-1= (I+PG)^-1 = S0

S0 = feedback(eye(2),F0);

%abbiamo poli a parte re<0
min(pole(S0))

%valori singolari
figure
sigma(S0,omega,'r')
hold on
grid on
sigma(inv(S0),omega,'b');
ps=sigma(inv(S0),omega);
legend('S0','ps');

ps=ps(1,:);
ps=gamma1*ps;
w1=make_approx_tf(ps,9,omega,'approx of ps in w1 scalar tf');


%% Calculate V0

%NB per vederli ad occhio basta sapere che V0 = -Wun ed T0 = -Wyn , quindi
%nello schema generale spegni i segnali che non ti servono e li becchi
%subito


close all;
%Semplificando lo schema esce una cosa del genere : 

%
%
%        (-)
%    n--->o--->G------->u
%         |         |
%         |____P<___|


%   
%                   
%
%

%Voglio che max val sing di V0<1/la 

V0 = feedback(LTR_control,nominal_model);
la_test=la.^-1;
figure(1);
max_V0=sigma(V0,omega,'r');
max_V0=max_V0(1,:);
hold on
grid on
semilogx(omega,mag2db(max_V0),'r');
semilogx(omega,mag2db(la_test),'b');
legend('V0','1/la');
title('V0< 1/la');


% 
% figure(2);
% semilogx(omega,mag2db(la))
% grid on;
% hold on;
% semilogx(omega,mag2db(la_signed))
% legend('la','la signed');
% title('la>la signed');

%calcolo approx w2 a partire da la
w2=make_approx_tf(la,9,omega,'approx of la in w2 scalar tf');

%w2 >=la
figure(2);
sigma(w2,omega,'r');
hold on;
grid on;
semilogx(omega,mag2db(la),'b')
legend('w2','la ');
title('w2>la');


%% Calcolo T0
close all;
T0 = feedback(F0,eye(2));
lm_test=lm.^-1;

figure(1);
max_T0=sigma(T0,omega);
max_T0=max_T0(1,:);
hold on
grid on
semilogx(omega,mag2db(max_T0),'r');
semilogx(omega,mag2db(lm_test),'b');
legend('T0','1/lm out')
title('T0 < 1/lm out');


%calcolo apporx w2 a partire da la
w3=make_approx_tf(lm,9,omega,'approx of lm in w3 scalar tf');

%w3 >=lm
figure(2);
sigma(w3,omega,'r');
hold on;
grid on;
semilogx(omega,mag2db(lm),'b')
legend('w3','1/lm out ');
title('w3>lm');

%% Check of pole and zeros of w1 w2 w3
disp('Zeri w1');
zero(w1)
disp('Zeri w2');
zero(w2)
disp('Zeri w3');
zero(w3)

disp('Poli w1');
pole(w1)
disp('Poli w2');
pole(w2)
disp('Poli w3');
pole(w3)


%% Synthesis of  H infinity 

epsilon=0;
plant = ss(An+epsilon*eye(length(An)),Bn,Cn,Dn)

eig(plant.A)

%% calculate expanded system

GAM_prev=3.1;
W1 =(1/GAM_prev)*w1*eye(2);
W2 =(1/GAM_prev)* w2*eye(2);
W3 = (1/GAM_prev)*w3*eye(2);

%la funzione augw ritorna il sistema allargato con le uscite fittizie w1 w2
%w3 viste a teoria per dare delle specifiche di robustezza e stabilità
P = augw(plant,W1,W2,W3);

%% calcolo (selezione) matrici sistema esteso 

[rB,cB] = size(P.B);
B1 = P.B(:,1:cB-2);
B2 = P.B(:,cB-1:cB);

[rC,cC] = size(P.C);
C1 = P.C(1:rC-2,:);
C2 = P.C(rC-1:rC,:);

D11 = P.D(1:rC-2,1:cB-2);
D12 = P.D(1:rC-2,cB-1:cB);
D21 = P.D(rC-1:rC,1:cB-2);
D22 = P.D(rC-1:rC,cB-1:cB);



%% Check condition of H infinity

%D22 = 0
D22

%D11 = 0 --> W1 strictly proper 
D11
figure;
sigma(W1,'o-');
title('Singular values of W1');
grid on;

%rg(D12) pieno colonna
D12
rank(D12)

%rg(D21) pieno riga
D21
rank(D21)

%% no null real part zeros of [A-sI,B2; C1, D12])

tzero(ss(P.A,B2,C1,D12))

%% no null real part zeros of [A-sI,B1; C2, D21]
disp(tzero(ss(P.A,B1,C2,D21)))

%% synthesis of controller 
[K,CL,GAM] = hinfsyn(P);
GAM

% definizione della catena diretta
F0 = series(K,ss_model_unc.NominalValue);


%% costruzione di S0 
S0 = feedback(eye(2),F0);
min(pole(S0))

%valori singolari
figure
sigma(S0,'r')
hold on
grid on
sigma(inv(W1),'g')
legend('S0','1/W1')

%% costruzione di V0
V0 = feedback(K,ss_model_unc.NominalValue);

%valori singolari
figure
sigma(V0,'r')
hold on
grid on
sigma(inv(W2),'g')
legend('V0','1/W2')

%% costruzione di T0
T0 = feedback(F0,eye(2));

%valori singolari
figure
sigma(T0,'r')
hold on
grid on
sigma(inv(W3),'g')
legend('T0','1/W3')


%% plot step of nominal system
figure
step(T0,'r')
grid on 
hold on
step(S0,'g')
legend('step(T0)','step(S0)')

%% make perturbation 
for i=1:1:N
   
    perturbation_tf_Min{i}=make_approx_tf(max_sig_dMin(i,:),9,omega,'');
    F_Temp=series(F0,eye+perturbation_tf_Min{i});
    perturbation_systems{i}=feedback(F_Temp,eye(2));    
end

%%
close all;
t_final=7;
row_of_plots=round(N/3);
figure(2);
for i=1:1:N
     subplot(row_of_plots,3,i);  
     str=strcat('Step Response of Variation : ', num2str(i));
     step(perturbation_systems{i},t_final,'k');
     title(str);


end


