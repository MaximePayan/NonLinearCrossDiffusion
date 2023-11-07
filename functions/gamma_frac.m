function y = gamma_frac(x,P,Q,nu,trunc)
    if class(P)~="Polynomial"
        P = Polynomial(P);
    end
    if class(Q)~="Polynomial"
        Q = Polynomial(Q);
    end
    if nargin == 4
        y = mult_err(polynom_err(x,P,nu), ...
    inv_err(polynom_err(x,Q,nu),nu),nu);
    else
        y = mult_err(polynom_err(x,P,nu,trunc), ...
    inv_err(polynom_err(x,Q,nu,trunc),nu),nu,trunc);
    end
end