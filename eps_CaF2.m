function eps = eps_CaF2(w)
    lam = 1./(w/10000);
    eps = 1.33973+(0.69913*lam.^2)./(lam.^2-0.09374^2)+(0.11994*lam.^2)./(lam.^2-21.18^2)+(4.35181*lam.^2)./(lam.^2-38.46^2);
    % eps = eps';
end