%% Clear all memory , close all plot and clear the command
clear all;
clc;
close all;
%% Load the model 

load('model.mat');
A=nominal_model.A
B=nominal_model.B;
C=nominal_model.C;
D=zeros(2,2)
system = ss(A,B,C,D);



%% Define the matrix of LQR(First synthesis)
Q =C'*C;
rho1 = 10;
R = rho1*eye(2);
    
%% First LQR's Synthesis


%Matrice Kopt matrice dei guadagni ottimi ( risolve problema ottima di LQR)
%La LQR ha forma u=Kx con K soluzione problema LQR con K=R^-1 * B' * P , e
%P soluzione della DRE. 

%Allora creiamo il sistema H(openLoop_nominal) con matrice C = Kopt
%che poi facendo CL=feedback(H,eye(2)) otteniamo che la u viene messa
%uguale alla y (chiudo anelli con guadagno unitario, passo eye(2)) e allora
%otteniamo che xdot = Ax + By ( u = y) , ma y=Cx[D=0] (ovvero
%C=Kopt)-->otteniamo xdot = (A+BKopt)*x che sarà asinotaticamente stabile
%dalla teorema di asintotica stabilità delle LQR 
% 
%           ---->[H]----->y
%           |          |                    xdot=Ax+Bu
%           |<---------|                    y=Cx+Du=Cx=u ( chiuso feedback        
%              feedback = I                        D=0     con I)
% 
%                                           xdot=Ax+By=Ax+B*Cx=Ax+B*kopt*x
%                                               =(A+B*Kopt)*x


Kopt1 = lqr(A,B,Q,R);

disp('Eig value of closed loop  ');
eigValue_K1 = eig(A-B*Kopt1)



%% Define matrix of Kalman filter for the first LQR controller
clc;
%Calcolo guadagno del filtro di Kalman
sigma1=10;
W=sigma1*(B*B');
V=eye(2);

%% Kalman's Synthesis for the first LQR controller
L1= lqr(A',C',W,V)';

disp(L1);
eig_Value_L1 = eig(A-L1*C);
disp(eig_Value_L1);


%% First contoller:  Regulator composed by LQR + Kalman
clc;
Ac = A-B*Kopt1-L1*C;
Bc = L1;
Cc = Kopt1;
Dc = zeros(size(Kopt1,1),size(L1,2));


controller1 = ss(Ac,Bc,Cc,Dc);
closedLoop1 = feedback(series(controller1, system),eye(2)); 
[Acl,Bcl,Ccl,Dcl] = ssdata(closedLoop1);

disp('Eig values of first closed loop system(LQR + Kalman) ');
eig_Value_CL1 = eig(Acl);
disp(sort(eig_Value_CL1));


figure(1);
step(closedLoop1);
title('Step response of first regulator');
%% New Synthesis (hopefull faster than first one) Second controller
clc;

sigma2=100;
rho2=100;
R=rho2*eye(2);
W=sigma2*(B*B');
V=eye(2);


Kopt2 = lqr(A,B,Q,R);
L2= lqr(A',C',W,V)';


eig_Value_K2= eig(A-B*Kopt2);
eig_Value_L2 = eig(A-L2*C);

disp('Eig values of two synthesis');
disp('K1         K2         L1         L2:');
disp([sort(real(eigValue_K1),'descend') sort(real(eig_Value_K2),'descend') sort(real(eig_Value_L1),'descend') sort(real(eig_Value_L2),'descend')]);

Ac = A-B*Kopt2-L2*C;
Bc = L2;
Cc = Kopt2;
Dc = zeros(size(Kopt2,1),size(L2,2));

controller2 = ss(Ac,Bc,Cc,Dc);
closedLoop2 = feedback(series(controller2,system),eye(2)); 

figure(1)
step(closedLoop2);
title('Step response of second regulator');



%% Shift the matrix A with original weight of LQR and Kalman
clc;

disp('Shift the A matrix with alpha=10 with original weight for matrix of LQR and Kalman')
Q =C'*C;

R = rho1*eye(2);
W=sigma1*(B*B');
V=eye(2);

alpha=10;


Kopt3 = lqr(A+alpha*eye(size(A)),B,Q,R);
L3= lqr(A'+alpha*eye(size(A)),C',W,V)';


eig_Value_K3=eig(A-B*Kopt3);
eig_Value_L3=eig(A-L3*C);

disp('Real part of eig values of the three synthesis:');
disp('   K1        K2         K3        L1         L2         L3:');
disp([sort(real(eigValue_K1),'descend') sort(real(eig_Value_K2),'descend') sort(real(eig_Value_K3),'descend') sort(real(eig_Value_L1),'descend') sort(real(eig_Value_L2),'descend') sort(real(eig_Value_L3),'descend')]);
disp('As you can see, the slowly eig values are shifted a lot to the left of Im axes');

Ac = A-B*Kopt3-L3*C;
Bc = L3;
Cc = Kopt3;
Dc = zeros(size(Kopt3,1),size(L3,2));


controller3 = ss(Ac,Bc,Cc,Dc);
closedLoop3 = feedback(series(controller3, system),eye(2));

figure(1);
step(closedLoop3);
title('Step response of third Controller (Shifted matrix)');

%% Synthesis with KLQG

alfa_system=ss(A+alpha*eye(size(A)),B,C,D); 

QXU=blkdiag(Q,R);
QWV=blkdiag(W,V);

KLQG=lqg(alfa_system,QXU,QWV);
KLQG.A=KLQG.A-alpha*eye(size(KLQG.A)); %Re-shifted regulator
disp('Eig values of closed loop with regulator made with KLQG');
Wcl=feedback(KLQG,system,1); %KLQG uses postive feedback
eig(Wcl)
figure(1)
step(Wcl);
title('Step response of KLQG controller');