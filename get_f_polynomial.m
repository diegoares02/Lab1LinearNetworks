function [Fs,Fw] = get_f_polynomial(R)
K=@(n)(sqrt(1-(1/(n^2))));
leftscroll = @(p1,p2) [zeros(1,max(0,numel(p2) - numel(p1))),p1];
u=[1 -1/R(1)];
v=[0 K(R(1))];
w_prima_cuadrada=[1 0 -1];
w=[1 0];
aux_u=u;
aux_v=v;
    for m=2:numel(R)
        p_u=K(R(m))*(conv(w_prima_cuadrada,aux_v));
        p_v=conv(w,aux_v);
        u=leftscroll(conv(w,aux_u),p_u)-leftscroll(aux_u/R(m),p_u)+p_u;
        v=p_v-leftscroll(aux_v/R(m),p_v)+leftscroll(K(R(m))*aux_u,p_v);
        aux_u=u;
        aux_v=v;
    end
Fw=poly(roots(u));
Fs=poly(roots(u)*1i);
display('F(w):');
display(Fw);
end