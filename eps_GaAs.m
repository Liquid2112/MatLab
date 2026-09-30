function eps = eps_GaAs(w) % GaAs
    wlGaAs      = 292;
    wtGaAs     = 268;
    GtGaAs     = 0.8;
    GlGaAs     = 0.8;
    epsinfGaAs = 10.9;
    eps =  epsilon_1phonon(w,wtGaAs,wlGaAs,GtGaAs,GlGaAs,epsinfGaAs);
end

function eps = epsilon_1phonon(w,wT,wL,gT,gL,epsinf)
    eps  = epsinf*(w.^2-wL^2+1i*gL*w)./(w.^2-wT^2+1i*gT*w);
end