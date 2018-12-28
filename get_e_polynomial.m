function[Ew,Es,Ew_h]=get_e_polynomial(Pw,Fw,epsilon)
leftscroll = @(p1,p2) [zeros(1,max(0,numel(p2) - numel(p1))),p1];
Ew=leftscroll((Pw/epsilon),Fw) - 1i*Fw;
Er=roots(Ew);%roots of E(w)

%Check Hurwitz condition
for n=1:numel(Er)
    if real(Er(n))<=0
       Ew_h(n)=Er(n);
    end
    if real(Er(n))>0
       Ew_h(n)=-real((Er(n)))+1i*imag(Er(n));
    end
end
Es=poly(Er*1i);
display(Er);
display(roots(Ew_h))
end