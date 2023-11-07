function [gammav,gammapv,gammappv,gammapppv,errgammav,errgammapv,errgammappv,errgammapppv,neighbor_errgammav, neighbor_errgammapv, neighbor_errgammappv] = transform_gamma(v,nu,r,func_opt,which_gamma,N,N_out,K)
if which_gamma == "frac_2"
    which_gamma = "frac";
end
switch which_gamma
    case 'exp'
        alpha = func_opt{end};
        [gammav,errgammav] = exp_err(alpha,v,K,N,N_out,nu);
        gammapv = alpha*gammav;

        if nargout > 2
            gammappv = alpha*gammapv;
            gammapppv = alpha*gammappv;
            errgammapv = abs(alpha)*errgammav;
            errgammappv = abs(alpha)*errgammapv;
            errgammapppv = abs(alpha)*errgammappv;
            alpha = func_opt{end};
            G1 = abs(alpha)*func_opt{1}((norme_nu(v,nu)+r));
            G2 = abs(alpha)*G1;
            G3 = abs(alpha)*G2;
            neighbor_errgammav = errgammav + G1*r;
            neighbor_errgammapv = errgammapv + G2*r;
            neighbor_errgammappv = errgammappv + G3*r;
        end
    
    case 'frac'
        P=func_opt{end-1};Q=func_opt{end};
        d1 = (max([length(P) length(Q)])-1)*(N-1)+1;
        d1 = max([d1 N_out]);
        one1 = zeros(d1,1);
        one1(1)=1;
        [Pv,Qv]=vectors_PvQv(v,P,Q,d1);
        inv_Qv = convomat(Qv)\one1;
        
        temp1 = convo(Pv,inv_Qv,d1);
        gammav = temp1(1:N_out);
        
        Pp = polynomial_derivative(P);
        Qp = polynomial_derivative(Q);
        Pnum = conv(Pp,Q) - conv(P,Qp);
        Pden = conv(Q,Q);
        d2 = (max([length(Pnum) length(Pden)])-1)*(N-1)+1;
        [Pnumv,Pdenv]=vectors_PvQv(v,Pnum,Pden,d2);
        one2 = zeros(d2,1);
        one2(1)=1;
        inv_Pdenv = convomat(Pdenv)\one2;
        temp2 = convo(Pnumv,inv_Pdenv,d2);
        gammapv = temp2(1:N_out);

        if nargout > 2
            
            err1 = norme_nu(convo(Qv,inv_Qv)-one1,nu); %it must be smaller than 1
            errgammav =  err1/(1-err1)*norme_nu(Pv,nu)*norme_nu(inv_Qv,nu);
            errgammav = errgammav + norme_nu(temp1(N_out+1:end),nu);

            bound_Pv = g_frac(norme_nu(v,nu)+r,Pp,1);
            bound_Qv = g_frac(norme_nu(v,nu)+r,Qp,1);
            [~,errinv_Qv] = inv_err(Qv,bound_Qv,nu);
            [~,neighbor_errgammav] = mult_err(Pv,bound_Pv,inv_Qv,errinv_Qv,nu); 

            err2 = norme_nu(convo(Pdenv,inv_Pdenv)-one2,nu); %it must be smaller than 1
            errgammapv =  err2/(1-err2)*norme_nu(Pnumv,nu)*norme_nu(inv_Pdenv,nu);
            errgammapv = errgammapv + norme_nu(temp2(N_out+1:end),nu);

            bound_Pnumv = g_frac(norme_nu(v,nu)+r,polynomial_derivative(Pnum),1);
            bound_Pdenv = g_frac(norme_nu(v,nu)+r,polynomial_derivative(Pden),1);
            [~,errinv_Pdenv] = inv_err(Pdenv,bound_Pdenv,nu);
            [~,neighbor_errgammapv] = mult_err(Pnumv,bound_Pnumv,inv_Pdenv,errinv_Pdenv,nu); 

            Ppp = polynomial_derivative(Pp);
            Qpp = polynomial_derivative(Qp);
            Q2 = conv(Q,Q);
            Qp2 = conv(Qp,Qp);
            Pnum1 = conv(Ppp,Q2); 
            Pnum21 = 2*conv(Pp,Qp); Pnum22 = conv(P,Qpp);
            n2 = max([length(Pnum21) length(Pnum22)]);
            Pnum2 = -conv(Q,[Pnum21,zeros(n2-length(Pnum21),1)]+[Pnum22,zeros(n2-length(Pnum22),1)]);
            Pnum3 = 2*conv(P,Qp2);
            n = max([length(Pnum1) length(Pnum2) length(Pnum3)]);
            Pnum = [Pnum1,zeros(1,n-length(Pnum1))]+[Pnum2,zeros(1,n-length(Pnum2))]...
                +[Pnum3,zeros(1,n-length(Pnum3))];
            Pden = conv(Q2,Q);
            d3 = (max([length(Pnum) length(Pden)])-1)*(N-1)+1;

            [Pnumv,Pdenv]=vectors_PvQv(v,Pnum,Pden,d3);
            [inv_Pdenv,errinv_Pdenv] = inv_err(Pdenv,0,nu);
            [gammappv,errgammappv] = mult_err(Pnumv,0,inv_Pdenv,errinv_Pdenv,nu,N_out);
            
            bound_Pnumv = g_frac(norme_nu(v,nu)+r,polynomial_derivative(Pnum),1);
            bound_Pdenv = g_frac(norme_nu(v,nu)+r,polynomial_derivative(Pden),1);
            [~,errinv_Pdenv] = inv_err(Pdenv,bound_Pdenv,nu);
            [~,neighbor_errgammappv] = mult_err(Pnumv,bound_Pnumv,inv_Pdenv,errinv_Pdenv,nu); 

            Pppp = polynomial_derivative(Ppp);
            Qppp = polynomial_derivative(Qpp);
            Q3 = conv(Q2,Q);
            Qp3 = conv(Qp2,Qp);
            Pnum1 = conv(Pppp,Q3);
            Pnum2 = 6*conv(conv(Q,Qp),conv(Pp,Qp)+conv(P,Qpp));
            Pnum31 = 3*conv(Ppp,Qp);
            Pnum32 = 3*conv(Pp,Qpp);
            Pnum33 = conv(P,Qppp);
            n3 = max([length(Pnum31) length(Pnum32) length(Pnum33)]);
            Pnum3 = -conv(Q2,[Pnum31,zeros(1,n3-length(Pnum31))]+[Pnum32,zeros(1,n3-length(Pnum32))]...
                +[Pnum33,zeros(1,n3-length(Pnum33))]);
            Pnum4 = -6*conv(P,Qp3);
            n = max([length(Pnum1) length(Pnum2) length(Pnum3) length(Pnum4)]);
            Pnum = [Pnum1,zeros(1,n-length(Pnum1))]+[Pnum2,zeros(1,n-length(Pnum2))]...
            +[Pnum3,zeros(1,n-length(Pnum3))]+[Pnum4,zeros(1,n-length(Pnum4))];
            Pden = conv(Q3,Q);

            d4 = (max([length(Pnum) length(Pden)])-1)*(N-1)+1;
            [Pnumv,Pdenv]=vectors_PvQv(v,Pnum,Pden,d4);
            [inv_Pdenv,errinv_Pdenv] = inv_err(Pdenv,0,nu);
            [gammapppv,errgammapppv] = mult_err(Pnumv,0,inv_Pdenv,errinv_Pdenv,nu,N_out);

        end
    case 'WX'
        [gammav,gammapv,gammappv,gammapppv,errgammav,errgammapv,errgammappv,errgammapppv,neighbor_errgammav, neighbor_errgammapv, neighbor_errgammappv] = transform_gWX(v,nu,r,N,N_out,K);
    otherwise
      disp('case not yet handled') 
end
end