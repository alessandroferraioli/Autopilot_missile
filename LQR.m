%% Clear all memory , close all plot and clear the command
clear all;
close all;
clc;
%% Load model and variations
load('model.mat');
load('variations.mat');
syms s;



%% Define nominal matrixs of system

An=nominal_model.A
Bn=nominal_model.B
Cn=nominal_model.C
Dn=nominal_model.D



%% Check for LQR 
close all;
clc;
disp('Design of the LQR controller');

%Define some matrixs
Q =Cn'*Cn;
P0=tf(nominal_model);

%Zeros check
zeros_P0=tzero(P0);
check=0;
for i =1 :1:length(zeros_P0)
   if(real(zeros_P0(i))==0)
       disp('Zero(s) exists on imaginary axes ');
       check=1;
   end
end

if(check==0)
    disp('Zero doesn t exist on imaginary axes ');
end

%% Define matrix weight of LQR
rho = 10;
R = rho*eye(2);

%% LQR's Synthesis
clc;
close all;

%Matrice Kopt matrice dei guadagni ottimi ( risolve problema ottima di LQR)
%La LQR ha forma u=Kx con K soluzione problema LQR con K=R^-1 * B' * P , e
%P soluzione della DRE. 

%Allora creiamo il sistema (openLoop_nominal) con matrice C = Kopt
%che poi facendo CL=feedback(openLoop_nominal,eye(2)) otteniamo che la u viene messa
%uguale alla y (chiudo anelli con guadagno unitario, passo eye(2)) e allora
%otteniamo che xdot = Ax + By ( u = y) , ma y=Cx[D=0] (ovvero
%C=Kopt)-->otteniamo xdot = (A+BKopt)*x che sarà asinotaticamente stabile
%dalla teorema di asintotica stabilità delle LQR 
% 
%           ---->[openLoop_nominal]----->y
%           |                      |        xdot=Ax+Bu
%           |<---------------------|        y=Cx+Du=Cx=u ( chiuso feedback        
%              feedback = I                 D=0     con I)
% 
%                                           xdot=Ax+By=Ax+B*Cx=Ax+B*kopt*x
%                                               =(A+B*Kopt)*x

Kopt = lqr(An,Bn,Q,R);

openLoop_nominal = tf(ss(An,Bn,Kopt,zeros(2,2)));
%Eigenvalues of closed loop have real part <0 
ss_nominal_closedLoop=feedback(test,eye(2));
disp('Eigenvalues of closed loop nominal model');
%eig(ss_nominal_closedLoop.A)

%Calculate closed loop transfer function
tf_nominal_closedLoop_LQR = (ss_nominal_closedLoop);
figure(1)
hold on
sigma(openLoop_nominal,'b')
sigma(tf_nominal_closedLoop_LQR,'k')
grid on
legend('OpenLoop','closedLoop')

figure(2);
grid on;
step(tf_nominal_closedLoop_LQR);
title('Step response of closedLoop ');

%% LTR 

%Kalman observator
clc;
disp('LTR recovery')

%sigma molto alto mi permette di ottenere prima la convergenza tra la
%funzione anello aperta nel caso in cui posso misurare lo stato e quella
%con il filtro di kalman


sigmas1=100;
sigmas2=5000;           


W1=sigmas1*(Bn*Bn');
W2=sigmas2*(Bn*Bn');

V=eye(2);

L1= lqr(An',Cn',W1,V)';
L2= lqr(An',Cn',W2,V)';

Ac1 = An-Bn*Kopt-L1*Cn;
Ac2 = An-Bn*Kopt-L2*Cn;

Bc1 = L1;
Bc2 = L2;

Cc = Kopt;
Dc = zeros(2);

control_1 = ss(Ac1,Bc1,Cc,Dc);
control_2 = ss(Ac2,Bc2,Cc,Dc);

ss_nominal_openLoop_LTR_1 = series(nominal_model,control_1);
ss_nominal_openLoop_LTR_2 = series(nominal_model,control_2);

sigma(ss_nominal_openLoop_LTR_1,'b')
hold on;
grid on;
sigma(ss_nominal_openLoop_LTR_2,'k')
sigma(openLoop_nominal,'r');
legend('LTR: sigma = 100',' LTR: sigma =3000',' LQR');

%% Calculate  upper bound lma_in(w)
clc;
close all;


tf_nominal_closedLoop_LTR=tf(feedback(ss_nominal_openLoop_LTR_2,eye(2)));

%U0_LTR approssima bene la U0
figure(1);
hold on;
sigma(tf_nominal_closedLoop_LQR,omega,'r');
grid on;
sigma(tf_nominal_closedLoop_LTR,omega,'b');
legend('LQR closedLoop','LTR closedLoop');
title(' LTR approximates really well LQR');

%make bound lma_signed(Basically i'm gonna use just that) 
temp=sigma(tf_nominal_closedLoop_LTR,omega);
lma_signed=temp(1,:);
lma_signed=lma_signed.^-1;

%make bound lma_signed_nominal, never used
temp=sigma(tf_nominal_closedLoop_LQR,omega);
lma_signed_nominal=temp(1,:);
lma_singed_nominal=lma_signed_nominal.^-1;

%----------------------------------------------------------------------------------
% %% Test difference between sigma(U0_LTR,omega,1) and (sigma(U0_LTR,omega))^-1
% close all;
% figure(1);
% sigma(tf_nominal_closedLoop_LTR,omega,1,'r');
% hold on;
% grid on;
% semilogx(omega,mag2db(lma_signed),'b');
% legend(' LTR type 1',' LTR ^ -1');
% disp(' Note the difference between type 1 and other');
% disp('In type 1 case , you re getting the minimun val sing (see theroy)');
%----------------------------------------------------------------------------------



%% lm_in<=lma_signed
%voglio che tutte le mie variazioni moltiplicative riportate sull'ingresso
%siano con max val sing <= lma_signed(il mio bound)
close all;
figure(1);
hold on;
semilogx(omega,mag2db(lma_signed),'b','LineWidth',2)%Bound
grid on;
Legend=cell(N,1);
Legend{1}='bound=1/max-val-sing(U0)';
for iter=2:N+1
   Legend{iter}=strcat('Variation', num2str(iter-1));
end
for i=1:1:N
 
  semilogx(omega,mag2db(max_sig_dMin(i,:)),'Color',rand(1,3),'LineWidth',1.5)
  temp_bound = frd(max_sig_dMin(i,:),omega);
  
  %Mi salvo le corrispondenti funzioni di trasferimento
  perturbation{i} = tf(fitmagfrd(temp_bound,9,[],[],1)); 

end

title('Max sing values: input multiplicative uncertainties (random), bound (blue)');
legend(Legend);

%Get variations out of bound
check=0;
for j=1:1:N %for each variations
        for i=1:1:length(omega) %for each frequency

        value_var=max_sig_dMin(j,i);%variazione j esima alla frequenza i
        value_bound=lma_signed(i);%bound frequenza i
        if(value_var>=value_bound)
            str_disp=sprintf( 'Variation %d is out of bound for frequency %f of %f db',j,i,mag2db(value_var-value_bound));
            disp(str_disp);
            check=1;
            break;
        end
        end
end
if(check==0)
    disp('All variations are below the bound');
end




%Ora dentro il mio oggetto perturbation ho le funzioni di trasferimento
%delle perturbazioni moltiplicative sull ingresso che stanno sotto il mio
%bound lm(calcolato come la funzione che meglio "fitta" il mio massim valor
%singolare U0^-1 contenuto nella variabile la

%------------------------------------------------------------------------------------------------------------------
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%------------------------------------------------------------------------------------------------------------------


%RICAPITOLANDO COSA E' STATO FATTO FINO AD ORA  :

% -Ho calcolato la LQR e la LTR per l'impianto nominale
% -Ho calcolato il mio upper bound di perturbazioni 'lm'(definito come
% max_val_sing(U0) ^-1
% -Ho calcolato la sua funzione di trasferimento di approx di ordine 9 chiamata lm
% -Ho fatto variare p della matrice A in un range--> mi ha dato delle variazioni sul processo
% -Ho preso N=9 di sample di variaizioni del sistema 
% -Ho calcolato le variazioni additive e da quelle le moltiplicative sull ingresso
% -Ho visto che le varaizioni stanno tutte sotto il mio bound
% -Avevo delle associazioni modulo-frequenza-->
%  ho calcolato le corrispondenti funzion di trasferimento e le ho messe nell'oggetto perturbation


%------------------------------------------------------------------------------------------------------------------
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%------------------------------------------------------------------------------------------------------------------


%% Robustness LTR sigma=100

close all;
disp('LTR recovery: robustness, sigma=100;')
final_time=5;
nominal_LTR_H1=feedback(ss_nominal_openLoop_LTR_1,eye(2));
for i=1:1:N
    perturbation_systems{i}=feedback(ss_nominal_openLoop_LTR_1,eye(2)+perturbation{i}*eye(2));
end

figure(1)
step(nominal_LTR_H1,'r--');
title('Step response of nominal model');
grid on;
row_of_plots=round(N/3);
figure(2);
for i=1:1:N
     subplot(row_of_plots,3,i);  
     str=strcat('Step Response with Variation : ', num2str(i));
     step(perturbation_systems{i},'k',5);
     title(str);


end


%% Robustezza LTR sigma 3000
disp('LTR recovery: robustness, sigma=3000;')

nominal_LTR_H2=feedback(ss_nominal_openLoop_LTR_2,eye(2));
for i=1:1:N
    perturbation_systems{i}=feedback(ss_nominal_openLoop_LTR_2,eye(2)+perturbation{i}*eye(2));
  
end

figure(1);
step(nominal_LTR_H2,'r--');
grid on;
title('Step response of nominal model');

figure(2);
for i=1:1:N
     subplot(row_of_plots,3,i);   
     str=strcat('Step Response with Variation : ', num2str(i));
     step(perturbation_systems{i},'k',final_time);
     title(str);

end


%% salvo dati per Hinfinit

save('lqr_ltr.mat','control_2');
