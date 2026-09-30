function eps = eps_Ge_THz(w) % epsilon V2O5 in model of 4-parameter oscillator
    % data from Gonzalo             
    eps = eps_Ge_polynomial(w);  
end

% function eps = epsilon_1phonon(w, wT, wL, gT, gL)
%     eps  = (wL.^2-w.^2-1i*gL*w)./(wT.^2-w.^2-1i*gT*w);
% end
% 
function eps = eps_Ge_polynomial(w) % epsilon SiO2 in model of 4-parameter oscillator
    l = 1./w.*10^4;
    n = 4 + 0.001106337 - 0.00314503./l + 0.492812./l.^2 - 0.601906./l.^3 + 0.982897./l.^4;
    eps = n.^2;
end