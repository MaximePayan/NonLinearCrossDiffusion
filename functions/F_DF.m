function [F,DF] = F_DF(u,v,gv,gpv,N,D,a,b,sigma,eps)

l = b-a;
if exist('intval','file')
    Lap = intval(laplacien(D,l));
    I = intval(eye(D));
    eps = intval(eps);
    sigma = intval(sigma);
    
else
    Lap = laplacien(D,l);
    I = eye(D);
end


Mu = convomat(u,D); %size DxD
Mgv = convomat(gv.App,D);
Mgpv = convomat(gpv.App,D);
uu = Mu*[u ; zeros(D-N,1)]; 
ugv = Mu*[gv.App;zeros(D-N,1)];


F1 = Lap.*ugv + sigma*([u;zeros(D-N,1)]-uu); %size D
F2 = eps*Lap.*[v;zeros(D-N,1)] + [u-v;zeros(D-N,1)]; %size D

F = [F1;F2]; %size 2D

DF = {diag(Lap)*Mgv + sigma*(I-2*Mu), diag(Lap)*Mgpv*Mu;
    I, eps*diag(Lap)-I}; %size 2Dx2D

end
