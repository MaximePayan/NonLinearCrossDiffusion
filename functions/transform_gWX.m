function [gammav,gammapv,gammappv,gammapppv,errgammav,errgammapv,errgammappv,errgammapppv, neighbor_errgammav, neighbor_errgammapv, neighbor_errgammappv] = transform_gWX(v,nu,r,N,N_out,K)
    one = 0*v;
    one(1) = 1;
    [gv,errgv] = exp_err(9,v-one,K,N,N_out,nu);
    gpv = 9*gv;
    one = 0*gv;
    one(1)=1;
    gv = one + gv;
    
    [gammav,errgammav] = inv_err(gv,errgv,nu);
    [g2v,errg2v] = mult_err(gv,errgv,gv,errgv,nu);
    [invg2v,errinvg2v] = inv_err(g2v,errg2v,nu);
    errgpv = 9*errgv;
    [gammapv,errgammapv] = mult_err(-gpv,errgpv,invg2v,errinvg2v,nu);
    if nargout>2
        gppv = 9*gpv;
        errgppv = 9*errgpv;
        gpppv = 9*gppv;
        errgpppv = 9*errgppv;
 
        [g3v,errg3v] = mult_err(g2v,errg2v,gv,errgv,nu);
        [invg3v,errinvg3v] = inv_err(g3v,errg3v,nu);
        [g2pv, errg2pv] = mult_err(gpv,errgpv,gpv,errgpv,nu);
        [gvgppv,errgvgppv] = mult_err(gv,errgv,gppv,errgppv,nu);
        numerator_gammappv = 2*g2pv-gvgppv;
        errnumerator_gammappv = 2*errg2pv + errgvgppv;
        [gammappv,errgammappv] = mult_err(numerator_gammappv,errnumerator_gammappv, ...
            invg3v,errinvg3v,nu);

        [g4v,errg4v] = mult_err(g3v,errg3v,gv,errgv,nu);
        [invg4v,errinvg4v] = inv_err(g4v,errg4v,nu);
        [g2vgpppv,errg2vgpppv] = mult_err(g2v,errg2v,gpppv,errgpppv,nu);
        [g3pv,errg3pv] = mult_err(g2pv,errg2pv,gpv,errgpv,nu);
        [gvgpvgppv,errgvgpvgppv] = mult_err(gvgppv,errgvgppv,gpv,errgpv,nu);
        numerator_gammapppv = -g2vgpppv - 6*g3pv + 6*gvgpvgppv;
        errnumerator_gammapppv = errg2vgpppv + 6*errg3pv + 6*errgvgpvgppv;
        [gammapppv,errgammapppv] = mult_err(numerator_gammapppv,errnumerator_gammapppv, ...
            invg4v,errinvg4v,nu);
        
        norm_v1 = norme_nu(v-one,nu);
        bound_gv = 9*exp(9*(norm_v1+r))*r;
        [~,neighbor_errgammav] = inv_err(gv,errgv+bound_gv,nu);

        bound_gpv = 9*bound_gv;
        bound_g2v = 18*(1+exp(9*(norm_v1+r)))*exp(9*(norm_v1+r))*r;
        [~,neighbor_errinvg2v] = inv_err(g2v,errg2v+bound_g2v,nu);
        [~,neighbor_errgammapv] = mult_err(gpv,errgpv+bound_gpv,invg2v,neighbor_errinvg2v,nu);
        
        bound_numerator_gammappv = 729*exp(9*(norm_v1+r))*(1+2*exp(9*(norm_v1+r)))*r;
        bound_g3v = 27*(1+exp(9*(norm_v1+r)))^2*exp(9*(norm_v1+r))*r;
        [~,neighbor_errinvg3v] = inv_err(g3v,errg3v+bound_g3v,nu);
        [~,neighbor_errgammappv] = mult_err(numerator_gammappv,errnumerator_gammappv+bound_numerator_gammappv, ...
            invg3v,neighbor_errinvg3v,nu);
    end
end
    