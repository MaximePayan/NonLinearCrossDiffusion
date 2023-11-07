%% Initialisation of the data
addpath functions functions/Polynomial functions/xAppErr
a = 0; b=4*pi; %boundary of the segment
eps = 1; %diffusion rate of the species
sigma = 0.6; %logistic parameter
l = b-a;
N = 50; %size of vectors for the Newton's Method
N_out = 100; %size of vectors for the validation
D = 6*N_out; %total size needed in in the build of the A matrix
r = 1e-6; %default radius to validate the numerical solution
nu = 1.0001; %parameter of the norm

%% Specific data for the case: gamma = 1/(1+exp(9(x-1))
alpha = 9;
K = 30;
WX_data = {alpha,K};
%% A start to get an approximate solution
if exist("Data/U_initial_WX.mat",'file')
    load("Data/U_initial_WX.mat")
else
    g1 = 1;
    gp1 = -2.25;
    perturb = 1e-3;
    U = guess_approx(a,b,sigma,eps,g1,gp1,N,perturb);
    clear g1 gp1 perturb
end
%% Newton's Method to get an approximate solution
it_max = 20;
tol = 1e-12;
which_gamma = 'WX';
FDF = @(X) F_DF_diffusion_croise(X,a,b,sigma,eps,nu,WX_data,which_gamma,N_out);
[U,E] = func_Newton([U(1:N);zeros(N_out-N,1);U(N+1:end);zeros(N_out-N,1)],FDF,it_max,tol);
%% First Computations
u = U(1:N_out);
v = U(N_out+1:end);

if exist('intval','file')
    u = intval(u); v = intval(v);
end

one = xAppErr(1);
fv = exp_err(alpha,v-eye(N_out,1),K,nu,N_out);
fv_r = exp_err(alpha,xAppErr(v,r)-one,K,nu,N_out);
gv = inv_err(one+fv,nu);
gv_r = inv_err(one+fv_r,nu);

gpv = mult_err(-alpha.*fv,power_err(gv,2,nu,N_out),nu,N_out);
gv_r2 = power_err(gv_r,2,nu,N_out);
gpv_r = mult_err(-alpha.*fv_r,gv_r2,nu,N_out);

gppv_r = alpha^2.*mult_err(mult_err(fv_r,fv_r-one,nu,N_out),...
    mult_err(gv_r2,gv_r,nu,N_out),nu,N_out);


[F,DF] = F_DF(u,v,gv,gpv,N_out,D,a,b,sigma,eps);

%% Proof of existence theorem

if exist('intval','file')
    fprintf('\nRigorous validation with interval arithmetic\n')
else
    fprintf('\nPrevalidation without interval arithmetic\n')
end

[Y,Z1,Z2,r_min,r_max] = proof(a,b,sigma,eps,nu,u,gv,gpv,r,gpv_r,gppv_r,F,DF);

%% Figure
fig = false;
if fig
    plot_cos(U(1:N_out),a,b,100,'r--',2)
    hold on
    plot_cos(U(N_out+1:2*N_out),a,b,100,'b',2)
    xlim([a,b])
    ylabel('$y$','Interpreter','latex','Rotation',0)
    xlabel('$x$', 'Interpreter', 'latex')
    legend('u','v')
    title("Solution")
    set(gca,'FontSize',20)
    axis square
end