function [F,DF] = F_DF_diffusion_croise(U,a,b,sigma,eps,nu,func_opt,which_gamma,N_out)

if class(U) == "xAppErr"
    U = U.App;
end

N = length(U)/2;
u = [U(1:N);zeros(N_out-N,1)];
v = [U(N+1:2*N);zeros(N_out-N,1)];
l = b-a;

Mu = convomat(u,N_out);
switch which_gamma
    case 'frac'
        gv = gamma_frac(v,func_opt{1},func_opt{2},nu,N_out);
        gpv = gamma_frac(v,func_opt{3},func_opt{4},nu,N_out);
    case 'WX'
        alpha = func_opt{1};
        K = func_opt{2};
        fv = exp_err(alpha,v-eye(N_out,1),K,nu,N_out);
        gv = inv_err(xAppErr(1)+fv,nu);
        gpv = mult_err(-alpha.*fv,power_err(gv,2,nu,N_out),nu,N_out);
    otherwise
        disp('The case has not yet been processed')
end
Mgv = convomat(gv.App);
Mgpv = convomat(gpv.App);
uu = Mu*u;
ugv = Mu*gv.App;
Lap = laplacien(N_out,l);
I = eye(N_out);

F1 = Lap.*ugv + sigma.*(u-uu);
F2 = eps.*Lap.*v +u-v;

F = [F1;F2];

DF = [diag(Lap)*Mgv + sigma*(I-2*Mu), diag(Lap)*Mgpv*Mu;
    I, eps*diag(Lap)-I];
end
