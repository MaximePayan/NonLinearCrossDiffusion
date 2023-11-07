function [Upp,E] = func_Newton(U,F_DF,it_max,tol)
    itt = 0;
    Upp=U;
    [FU,DFU] = F_DF(Upp);
    e=norm(FU,1);
    E=e;
    while e>tol && itt<it_max
        Upp = Upp - DFU\FU;
        [FU,DFU]=F_DF(Upp);
        e = norm(FU,1);
        E=[E,e];
        itt = itt + 1;
    end
end

