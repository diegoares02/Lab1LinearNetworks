K=@(n)(sqrt(1-(1/n^2)));
leftpadz = @(p1,p2) [zeros(1,max(0,numel(p2) - numel(p1))),p1];

R=[1.3217 1.8082];%W Transmission zero roots
%R=[-1.7 -1.4 1.4 1.7];
N=4;
RL=20;
%P Parameters%
pw=poly(R);
%for m=2:numel(R)
%    p=conv(p,[1 -1i*R(m)]);
%    pw=conv(pw,[1 -R(m)]);
%end
if mod((N-numel(R)),2)==0
        ps=poly(1i*R)*1i;
    else
        ps=poly(R*1i);
end
%F Parameters%
u=[1 -1/R(1)];
v=[0 K(R(1))];
%w_prima_square=[1 0 -1];
aux_u=u;
aux_v=v;
%w=[1 0];
for m=2:numel(R)
    p_u=K(R(m))*(conv([1 0 -1],aux_v));
    p_v=conv([1 0],aux_v);
    u=leftpadz(conv([1 0],aux_u),p_u)-leftpadz(aux_u/R(m),p_u)+p_u;
    v=p_v-leftpadz(aux_v/R(m),p_v)+leftpadz(K(R(m))*aux_u,p_v);
    aux_u=u;
    aux_v=v;
end
for n=1:N-numel(R)
    p_u=(conv([1 0 -1],aux_v));
    p_v=conv([1 0],aux_v);
    u=leftpadz(conv([1 0],aux_u),p_u)+p_u;
    v=p_v+leftpadz(aux_u,p_v);
    aux_u=u;
    aux_v=v;
end
Fw=u;
Fs=poly(roots(Fw)*1i);
%Epsilon%
p_sum=0;
u_sum=0;
for m=1:numel(pw)
    p_sum=p_sum+pw(m);    
end
for m=1:numel(u)
    u_sum=u_sum+u(m);    
end
p_sum=abs(p_sum);
u_sum=abs(u_sum);
epsilon_r=1;
epsilon=sqrt(10^(RL/10))^(-1)*p_sum/u_sum;

%Calculo E(w)
E=leftpadz((pw/epsilon),u) - 1i*u;
Er=roots(E);%roots of E(w)
%for n=1:numel(Er)
%    if imag(Er(n))>0
%        Er(n)=conj(Er(n));
%    end
%end
Es=poly(Er*1i);
Es_root=roots(Es);
plot(Er,'x','Linewidth',2)
grid on
title('E(s)');
xlabel('Re');
ylabel('Im');



%S11 y S21

%d=(-5:0.001:5);
%s11=polyval(Fs,1i*d)./polyval(Es,1i*d);
%s21=polyval(ps,1i*d)./polyval(Es,1i*d)/epsilon;

%plot(d,20*log10(abs(s21)))
%hold on
%grid on
%plot(d,20*log10(abs(s11)),'r')
