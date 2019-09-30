%% Make an approximation of the sing_value from the frequency with specified order
function approx_tf = make_approx_tf(sing_value,order,frequency,txt)
pre_bound = frd(sing_value,frequency);
bound = fitmagfrd(pre_bound,order,[],[],1); 
bound_order = sigma(bound,frequency);

figure(100);
semilogx(frequency,mag2db(sing_value))
grid on;
hold on;
semilogx(frequency,mag2db(bound_order(1,:)));
legend('sing value','approx');
title(txt);

approx_tf=tf(bound);



end