function eps = eps_MgO_exp(w) % epsilon MgO in model of 4-parameter oscillator           
    eps = eps_MgO_3parOsc(w);  
end

function eps1 = epsilon1_phonon(w, S, wTO, gamma)
    eps1  = (S*wTO.^2*(wTO.^2-w.^2))./((wTO.^2-w.^2).^2+gamma.^2*w.^2);
end

function eps2 = epsilon2_phonon(w, S, wTO, gamma)
    eps2  = (S*wTO.^2*gamma*w)./((wTO.^2-w.^2).^2+gamma.^2*w.^2);
end

function eps = eps_MgO_3parOsc(w) % epsilon MgO in model of 4-parameter oscillator
    epsinf = 3.018;
    eps1 = zeros(size(w,1),size(w,2));
    eps2 = zeros(size(w,1),size(w,2));

    gamma = [18, 20];
    wTO = [403, 416]; 
    S = [1.5, 3.9];   

    for i = 1:length(wTO)
        eps1 = eps1 + epsilon1_phonon(w, S(i), wTO(i), gamma(i));
        eps2 = eps2 + epsilon2_phonon(w, S(i), wTO(i), gamma(i));
    end   
    eps = epsinf + eps1 + 1i*eps2;
end