leftpadz = @(p1,p2) [zeros(1,max(0,numel(p2) - numel(p1))),p1];

R=[-1.7 -1.4 1.4 1.7]; %Values of the finite trasmission zeros
N=8; %N=order of the filter
RL=20; %RL=return loss or reflection loss

%How to obtain the P(s) parameters
if isreal(R)
   p_s=poly(R)*1i;
else
    p_s=poly(R);
end

if mod((N-length(R)),2)==0
    p_s=poly(1i*R)*1i;
else 
    p_s=poly(1i*R);
end

raices_p_s=roots(p_s);
disp('The roots of P(s) are: ');
disp(raices_p_s);
disp('The coefficients of P(s) are: ');
disp(p_s);

%How to obtain the P(w) parameters
if isreal(R)
   p_w=poly(R);
else
    p_w=poly(R*-1i);
end

raices_p_w=roots(p_w);
disp('The roots of P(w) are: ');
disp(raices_p_w);
disp('The coefficients of P(w) are: ');
disp(p_w);

%How to obtain the F Parameters: 
%we can use the recursive method to calculate F(w) since we know that 
%F(w)=Cn(W) where n is the order of the filter
u1=[1 -1/R(1)];
v1=sqrt(1-1/R(1)^2);
w_prima_cuadrado=[1 0 -1];
w=[1 0];
aux_u=u1;
aux_v=v1;

%The for cycle when we have finite trasmission zeros
for m=2:length(R)
    k_u=conv(w_prima_cuadrado,aux_v)*sqrt(1-1/R(m)^2);
    k_v=conv(w,aux_v);
    y_v=sqrt(1-1/R(m)^2)*aux_u;
    u=leftpadz(conv(w,aux_u),k_u)-leftpadz(aux_u/R(m),k_u)+k_u;
    v=k_v-leftpadz(aux_v/R(m),k_v)+leftpadz(y_v,k_v);
    aux_v=v;
    aux_u=u;
end

%The for cycle for the infinte trasmission zeros
for n=1:N-length(R)
    k_u=(conv(w_prima_cuadrado,aux_v));
    k_v=conv(w,aux_v);
    u=leftpadz(conv(w,aux_u),k_u)+k_u;
    v=k_v+leftpadz(aux_u,k_v);
    aux_v=v;
    aux_u=u;
end
raices_f_w=roots(u);
f_w=poly(raices_f_w); %normalized f_w
f_s=poly(raices_f_w*1i);
disp('The coefficients of F(w) are: ');
disp(f_w);

%How to obtain Epsilon
if length(R)==N
       X=(1/sqrt((10^(RL/10))-1))*abs(sum(p_w)/sum(f_w));
       epsilon=sqrt(X^2+1);
       epsilon_r=epsilon/sqrt(epsilon^2-1);
else
       epsilon_r=1;
       epsilon=(1/(sqrt(10^(RL/10)-1)))*(sum(p_w)/sum(f_w)); 
end 
disp('The value of epsilon is:');
disp(epsilon);

%How to obtain E(w)
e_w =leftpadz(p_w/epsilon,f_w)-1i*f_w/epsilon_r;
raices_e_w= roots(e_w);

%We know that E(w) is a Hurwitz polynomial so all its roots must 
%lie on the left half part of the plane. So we need to convert the sign of
%all the roots fine before that present a negative Imaginary sign. In fact
%we are looking if the Imaginary part of the root is less than 0 and in
%case affermative the take the coniugate with the function conj.

 e_s=poly(raices_e_w*1i);
 raices_e_s=roots(e_s);
 
 figure(3)
 plot(raices_e_s,'rx');
 hold on; grid on;
 
for d=1:length(raices_e_s); 
    if real(raices_e_s(d))<=0
        raices_e_sH(d)=raices_e_s(d);
    end
    if real(raices_e_s(d))>0
        raices_e_sH(d)=-real(raices_e_s(d))+1i*imag(raices_e_s(d));
    end
end 

e_sH=poly(raices_e_sH);
disp('The roots of E(w) are:');
disp(raices_e_w);
disp('The coefficients of E(w) are: ');
disp(e_w);

d=(-5:0.001:5);
s11=polyval(f_s ,1i*d)./polyval(e_s, 1i*d);
s21=polyval(p_s ,1i*d)./polyval(e_s, 1i*d)/epsilon;

figure(1)
plot(raices_e_sH,'bx');
hold on;
grid on;

figure(2)
plot(d, 20*log10(abs(s21)),'blue'); hold on; grid on;
plot(d, 20*log10(abs(s11)),'red'); hold on; grid on;







