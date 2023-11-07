function U = guess_approx(a,b,sigma,eps,g1,gp1,N,perturb)
%inputs : parameters of the cross diffusion problem (a,b,sigma,eps,gamma,gammap),
%N the size of Fourier's vectors,
%perturb is a coeff of perturbation arround the stationnary state
%outputs : U=[u;v] of size 2*N which gives a start point for Newton's
%method, this point is 'near' the stationnary states (1,1), all this
%computation is based on the fact that (1,1) is a solution of the system.
l=b-a;
c1 = g1+gp1+eps*sigma;
c2 = 4*eps*sigma*g1;
c3 = (1-sigma)*(eps-g1);

K=[];
% K contains the k for which the Fourier's mode is amplified
if gp1*(c3-gp1)>=0 || c3-2*gp1>=0 
%     conditions for existence of a potential positive eigenvalue
    if c1^2 >= c2 && c1<0
%         condtion for those eigenvalues to be effectively positives
        b1 = (-c1 - sqrt(c1^2-c2))/(2*eps*g1);
        b2 = (-c1 + sqrt(c1^2-c2))/(2*eps*g1);
    else 
        b1=1;b2=-1;
        disp('eigenvalues are all negatives for all k')
    end
    for k=0:N-1
        test = (k*pi/l)^2;
        if test >= b1 && test <= b2
            K=[K,k];
        end
    end
else
    disp("all eigenvalues are complexes with negative real part or eigenvalues are negatives, for all k")
end
%d1 = 1+sigma +(g1+eps)*(K.*pi./l).^2;
delta = (eps-g1)^2*(K.*pi./l).^4 + 2*((sigma-1)*(g1-eps)-2*gp1)*(K.*pi./l).^2 + (sigma-1)^2 ;
%lamda_plus = 0.5*(-d1+sqrt(delta));
d2 = 1-sigma +(eps-g1)*(K.*pi./l).^2;
x = 0.5.*(d2 + sqrt(delta));
y = 1 + 0.*K;
u = zeros(N,1);
v = zeros(N,1);
u(K+1) = perturb.*x;
v(K+1) = perturb.*y;
u(1)=1;v(1)=1;
U=[u;v];
end

