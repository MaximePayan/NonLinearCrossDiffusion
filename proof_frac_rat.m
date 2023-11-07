%% Initialisation of the data
addpath functions functions/Polynomial functions/xAppErr
a = 0; b=3*pi; %boundary of the segment
eps = 1; %diffusion rate of the species
sigma = 0.053; %logistic parameter
l = b-a;
N = 50; %size of vectors for the Newton's Method
N_out = 100; %size of vectors for the validation
D = 6*N_out; %total size needed in in the build of the A matrix
r = 1e-6; %default radius to validate the numerical solution
nu = 1.0001; %parameter of the norm


%% Precalculations for the rational fraction case

P = Polynomial(1);
Q = Polynomial([1 0 0 0 0 0 0 0 0 1]);

Pp = derivative(P);
Qp = derivative(Q);
Pnum1 = Pp*Q - P*Qp;
Q2 = Q*Q;

Ppp = derivative(Pp);
Qpp = derivative(Qp);
Pnum2 = Q2*Ppp - Q*(2.*Pp*Qp+P*Qpp) + 2.*P*Qp*Qp;
Q3 = Q2*Q;

Pppp = derivative(Ppp);
Qppp = derivative(Qpp);

Pnum3 = Pppp*Q3 + 6.*Q*Qp*(Pp*Qp + P*Qpp) - ...
    Q2*(3.*Ppp*Qp + P*Qppp) - 6.*P*Qp*Qp*Qp;
Q4 = Q3*Q;

frac_data = {P,Q,Pnum1,Q2,Pnum2,Q3,Pnum3,Q4};
clear Pp Ppp Pppp Qp Qpp Qppp P Q Pnum1 Q2 Pnum2 Q3 Pnum3 Q4
%% A start to get an approximate solution
if exist("Data/U_initial_frac_rat.mat",'file')
    load("Data/U_initial_frac_rat.mat")
else
    g1 = scalar_eval(1,frac_data{1})/scalar_eval(1,frac_data{2});
    gp1 = scalar_eval(1,frac_data{3})/scalar_eval(1,frac_data{4});
    perturb = 1e-3;
    U = guess_approx(a,b,sigma,eps,g1,gp1,N,perturb);
    clear g1 gp1 perturb
end
%% Newton's Method to get an approximate solution

it_max = 20;
tol = 1e-12;
which_gamma = 'frac';
FDF = @(X) F_DF_diffusion_croise(X,a,b,sigma,eps,nu,frac_data,which_gamma,N_out);
U = func_Newton([U(1:N);zeros(N_out-N,1);U(N+1:end);zeros(N_out-N,1)],FDF,it_max,tol);
%% First Computations
u = U(1:N_out);
v = U(N_out+1:end);

if exist('intval','file')
    u = intval(u); v = intval(v);
end


gv = gamma_frac(v,frac_data{1},frac_data{2},nu,N_out);
gpv = gamma_frac(v,frac_data{3},frac_data{4},nu,N_out);
gpv_r = gamma_frac(xAppErr(v,r),frac_data{3},frac_data{4},nu,N_out);
gppv_r = gamma_frac(xAppErr(v,r),frac_data{5},frac_data{6},nu,N_out);

[F,DF] = F_DF(u,v,gv,gpv,N_out,D,a,b,sigma,eps);

%% Proof of existence theorem

if exist('intval','file')
    fprintf('\nRigorous validation with interval arithmetic\n')
else
    fprintf('\nPrevalidation without interval arithmetic\n')
end

[Y,Z1,Z2,r_min,r_max] = proof(a,b,sigma,eps,nu,u,gv,gpv,r,gpv_r,gppv_r,F,DF);