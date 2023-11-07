%% Initialisation of the data
addpath functions functions/Polynomial functions/xAppErr
a = 0; b=3*pi; %boundary of the segment
eps = 1; %diffusion rate of the species
l = b-a;
N = 50; %size of vectors for the Newton's Method
N_out = 150; %size of vectors for the validation
D = 6*N_out; %total size needed in in the build of the A matrix
r = 1e-8; %default radius to validate the numerical solution
nu = 1.0001; %parameter of the norm

load("Data/UU_discrete_diagram.mat")

Sigma = UU_refine(1,1:end);

%Initialisation of the list of index where the theorem
Itrue=[]; %concludes there exists solution
Ifalse=[]; %deos not conclude
Inocvg=[]; %says that the neighborhood is not good.

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
%% Calculations on each approximate solution in UU_refine
for ii=1:length(II)
    nbpoint = II(ii); %We apply the method on the current numerical solution UU_refine
    U_refine = UU_refine(2:end,nbpoint);
    sigma = Sigma(nbpoint);
    u = U_refine(1:N_out);
    v = U_refine(N_out+1:end);
    if exist('intval','file')
        u = intval(u); v = intval(v);
    end
    gv = gamma_frac(v,frac_data{1},frac_data{2},nu,N_out);
    gpv = gamma_frac(v,frac_data{3},frac_data{4},nu,N_out);
    gpv_r = gamma_frac(xAppErr(v,r),frac_data{3},frac_data{4},nu,N_out);
    gppv_r = gamma_frac(xAppErr(v,r),frac_data{5},frac_data{6},nu,N_out);
    [F,DF] = F_DF(u,v,gv,gpv,N_out,D,a,b,sigma,eps);
    [Y,Z1,Z2,r_min,r_max] = proof(a,b,sigma,eps,nu,u,gv,gpv,r,gpv_r,gppv_r,F,DF);
    %check hypothesis
    if not(Z1 < 1 && 2*Y*Z2 < (1-Z1)^2) %constants, are they good ?
        Ifalse = [Ifalse;nbpoint];
    elseif not(abs(r-(1-Z1)/Z2)<sqrt((1-Z1)^2-2*Y*Z2))%neighboorhood, is it good ? 
        Inocvg = [Inocvg;nbpoint];
    else %if the answers are yes ! The theorem is true for this one !
        Itrue = [Itrue;nbpoint];
    end
    disp(ii)
end
%% Figure
fig = false;
if fig
UU0_refine = UU_refine(2,:)+2*sum(UU_refine(3:N+1,:),1);
% Bifurcation is obtained in the numerical continuation, it contains the
% bifurcations numerically detected
B0 = Bifurcations(2,:) + 2*sum(Bifurcations(3:N+1,:),1);
Bs = Bifurcations(1,abs(B0)<10);
B0 = B0(abs(B0)<10);

% bifurc_points_11 gives the values of sigma for which (1,1) is a 
% bifurcation
g1 = scalar_eval(1,frac_data{1})/scalar_eval(1,frac_data{2});
gp1 = scalar_eval(1,frac_data{3})/scalar_eval(1,frac_data{4});
Spoints = bifurc_points_11(g1,gp1,eps,l,false);

figure(2)
clf(2,'reset')
plot(Sigma(Itrue),UU0_refine(Itrue),Color=[0.47,0.67,0.19],Marker='.', ...
    LineStyle='None',LineWidth=0.1,MarkerSize=9,DisplayName='validated solution')
hold on
plot(Bs,B0,Color=[0.49,0.18,0.56],Marker='o',LineStyle='None', ...
    LineWidth=2,MarkerSize=10, DisplayName='bifurcation numerically detected')
if exist('intval','file')
    for i=1:length(Spoints)
        plot([Spoints(i).inf,Spoints(i).sup],[1,1],Color=[1.00,0.41,0.16],Marker='None',...
            LineStyle='-',LineWidth=2,HandleVisibility='off')
    end
    plot(Spoints.mid,ones(size(Spoints)),Color=[1.00,0.41,0.16],Marker='*',...
        LineStyle='None',LineWidth=2,MarkerSize=10,DisplayName='bifurcation theoreticaly detected')
else
    plot(Spoints,ones(size(Spoints)),Color=[1.00,0.41,0.16],Marker='*',...
        LineStyle='None',LineWidth=2,MarkerSize=10,DisplayName='bifurcation theoreticaly detected')
end
legend
xlabel('\sigma','Rotation',0,'FontSize',18,'VerticalAlignment','bottom');
ylabel('u(0)','Rotation',0,'FontSize',18,'VerticalAlignment','cap');
axis([0 0.65 0.4 3.7])
end