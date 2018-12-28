function [Pw,Ps]=get_p_polynomial(R,N)
    if isreal(R)
        Pw=poly(R);
        Ps=poly(R)*1i;
    else
        Pw=poly(R*-1i);
        Ps=poly(R);
    end
    %check orthogonality
    if mod((N-numel(R)),2)==0
        Ps=poly(1i*R)*1i;
    else
        Ps=poly(R*1i);
    end
    display('P(w):');
    display(Pw);
end