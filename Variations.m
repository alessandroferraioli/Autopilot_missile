%% Clear all memory , close all plot and clear the command
clear all;
close all;
clc;

%% Load model
 
load('model.mat');
syms s;


%% Define nominal model
An=nominal_model.A
Bn=nominal_model.B
Cn=nominal_model.C
Dn=nominal_model.D

%% Make variations from model
clc;
close all;

% Get samples from uncertain model and make variations for each sample 
for i=1:1:N
%prendo campioni dal modello incerto 
sys{i} = usample(ss_model_unc);

%variazione additiva 
deltaA_sys{i} = tf(sys{i}) - tf(nominal_model);

%variazione moltiplicativa riportata sull' ingresso
deltaMin_sys{i} = (inv(tf(nominal_model))) * deltaA_sys{i};

%variazione moltiplicativa riportata sull' uscita
deltaMout_sys{i} = deltaA_sys{i} * (inv(tf(nominal_model)));

end



%get max singular values of perturbation
for i=1:1:N
  %sys{i} è un modello nello spazio di stato A B C D
  %sys  era il campione generato da usample a partire dal modello uss(A B C
  %D) con uss
  
  temp = sigma(sys{i},omega);
  max_sig_unc(i,:) = temp(1,:);
  
  temp = sigma(deltaA_sys{i},omega);
  max_sig_dA(i,:) = temp(1,:);
  
  temp = sigma(deltaMin_sys{i},omega);
  max_sig_dMin(i,:) = temp(1,:);
  
  temp = sigma(deltaMout_sys{i},omega);
  max_sig_dMout(i,:) = temp(1,:);
end

%define upper bound of variations
top_unc = max(max_sig_unc);
top_dA = max(max_sig_dA);
top_dMin = max(max_sig_dMin);
top_dMout = max(max_sig_dMout);

%save all
save('variations.mat','top_unc','top_dA','top_dMin','top_dMout','max_sig_unc','max_sig_dA','max_sig_dMin','max_sig_dMout');


