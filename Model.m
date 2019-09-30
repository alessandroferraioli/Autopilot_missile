%% Clear all memory , close all plot and clear the command

clear all;
close all;
clc;
%% Make model
%variabile p parametro messa a 17.5(valore nominale) come designato dalle schede
global A;
global B;
global C;
global D;
global ss_model_unc;
global nominal_model;
global p;
global omega;
global N;

%Number of samples
N=9;

%Frequency array
omega = logspace(-2,1,150);

%Range of variation for p parameter
p = ureal('p',17.5,'Range',[16.5 18.5]); 

%Declare the matrixs
A=[
   [3.23,p,-0.476,0,0.228,0,0,0];
   [-p,-3.23,0,0.476,0,-0.228,0,0];
   [0.39,0,-1.93,-0.8*p,-0.415,0,0,0];
   [0,-0.39,0.8*p,-93,0,-0.415,0,0];
   [0,0,0,0,0,0,0.75,0];
   [0,0,0,0,0,0,0,-0.75];
   [0,0,22.4,0,-0.300,0,-0.150,0];
   [0,0,0,-22.4,0,0.300,0,-0.150]
   ];
B=[[0,0,0,0,0,0,-1,0];[0,0,0,0,0,0,0,-1]]';
C=[[-2.99,0,-1.19,0.123*17.5,-27.64,0,0,0];[0,-2.99,0.123*17.5,1.19,0,27.64,0,0]];
D=0;

%make uncertain state space syste,
ss_model_unc=uss(A,B,C,D);

%Nominal  model
nominal_model=ss_model_unc.NominalValue;

%A has instable eigenvalues
disp('Eigenvalues of the matrix A')
disp(eig(A.NominalValue))

disp('Transfer matrix of the model')
tf_nominal_model=tf(nominal_model)

%As i said, the system results instable
figure(1);
step(tf_nominal_model);
grid on;
title("Step response of nominal model(without any kind of control)-->unstable");
% We won't  use it 
dim_system=size(A);

%Save model
save model.mat

%% Check the reachbility and observability of the system
disp('Check Reachability and Observability ')
disp('');

%Save nominal matrix(just to make more confortable the check)
nominalA=nominal_model.A;
nominalB=nominal_model.B;
nominalC=nominal_model.C;
ctrb(nominalA,nominalB);
if(rank(ctrb(nominalA,nominalB))==dim_system(1))
    disp('The system is reachable')
else 
    disp('The system isn t reachable')
 
end

if(rank(ctrb(nominalA',nominalC'))==dim_system(1))
    disp('The system is observable')
else 
    disp('The system isn t observable')
end

figure(1)
sigma(tf_nominal_model);
grid on
title('max and min singular value of nominal model');

