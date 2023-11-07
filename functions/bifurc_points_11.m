function Spoints=bifurc_points_11(f1,fp1,eps,l,fig)
%%% Return the sigmas for which a eigen value changes its sign, meansn it
%%% detectes bifurcation points around the state (1,1).

TrA = @(k,s) -(1+s+(f1+eps).*(k*pi/l).^2);
detA = @(k,s) eps*f1*(k*pi/l).^4 + (fp1+f1+eps*s).*(k*pi/l).^2 + s;
Delta = @(k,s) TrA(k,s).^2 - detA(k,s);

lambda = @(k,s) 0.5*(TrA(k,s)+sqrt(Delta(k,s)));

% K = 0:1:8;
% X = 0:0.1:8;
% sigma = 0.6;
% 
% Ypoints = lambda(K,sigma);
% Y = lambda(X,sigma);
% 
% figure(1)
% clf(1,'reset')
% plot(X,Y,'b-',K,Ypoints,'bo',K,0*K,'r-');
% title('\lambda_{+}');xlabel('k');ylabel('\lambda_{+}');

sigma_crit = @(k) -((eps*f1*(k*pi/l).^2 + (fp1+f1)).*(k*pi/l).^2)./(eps*(k*pi/l).^2+1);

k=1;
Spoints = [];
while sigma_crit(k)>=0
    Spoints = [Spoints,sigma_crit(k)];
    k=k+1;
end

if fig == true
    figure(2)
    clf(2,'reset')
    S = linspace(0,max(Spoints),100);
    X = linspace(0,length(Spoints)+1,100);
    [XX,SS]=meshgrid(X,S);
    Lambda_map = lambda(XX,SS);
    surf(XX,SS,Lambda_map)
    hold on
    surf(XX,SS,Lambda_map*0,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','cyan')
    for k=1:length(Spoints)
        sk = sigma_crit(k);
        plot3(k,sk,lambda(k,sk),'r',Marker="*",MarkerSize=15,LineWidth=3,DisplayName=string(k))
        text(k-0.02*abs(max(X)),sk,lambda(k,sk)+0.15*abs(max(Lambda_map,[],'all')),'('+string(k)+','+string(sk)+')')
        hold on
    end
    legend('\lambda^+','z=0','\sigma_k^b')
    title('\lambda^{+}_k');ylabel('\sigma');xlabel('k');zlabel('\lambda^{+}');
end
end
