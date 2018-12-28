clear;
clc;
leftscroll = @(p1,p2) [zeros(1,max(0,numel(p2) - numel(p1))),p1];
K=@(n)(sqrt(1-(1/n^2)));
% N=8;%Order of the filter
% RL=20;%Return Loss
% R=[-1.7 -1.4 1.4 1.7];%W Transmission zero roots
%-------------------------------------------------------------------------%
display('Lab Session 1');
%-------------------------------------------------------------------------%
display('1. Variables');
N=input('input the order of the filter N: \n');
RL=input('input the return loss of the filter RL: \n');
R=input('input the transmission zeros of the filter [1 2 3 4 5]: \n');
%-------------------------------------------------------------------------%
display('2. Calculation of polynomials');
NumInfRoots=N-numel(R);
R=[R inf(1,NumInfRoots)];
[Pw,Ps]=get_p_polynomial(R,N);
[Fs,Fw]=get_f_polynomial(R);
[epsilon,epsilon_r]=get_epsilon(R,N,Pw,Fw,RL);
[Ew,Es,Ew_h]=get_e_polynomial(Pw,Fw,epsilon);
%-------------------------------------------------------------------------%
display('3. Plot');
figure
subplot(3,3,[1 2 4 5 7 8]);
%S11 y S21
d=(-5:0.001:5);
s11=polyval(Fs,1i*d)./polyval(Es,1i*d);
s21=polyval(Ps,1i*d)./polyval(Es,1i*d)/epsilon;

plot(d,20*log10(abs(s21)))
hold on
grid on
plot(d,20*log10(abs(s11)),'r')
title('S-Parameters');

subplot(3,3,3)
plot(roots(Ew),'x','Linewidth',2)
grid on
title('E(w)');
xlabel('Re');
ylabel('Im');

subplot(3,3,9)
plot(Ew_h,'x','Linewidth',2)
grid on
title('E Hurwitz Condition');
xlabel('Re');
ylabel('Im');



