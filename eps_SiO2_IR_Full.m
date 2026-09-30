function eps = eps_SiO2_IR_Full(w) % epsilon V2O5 in model of 4-parameter oscillator
    % data from Gonzalo             
    eps = eps_SiO2_3parOsc(w);  
end

function eps = epsilon_1phonon(w, wT, wL, gT, gL)
    eps  = (wL.^2-w.^2-1i*gL*w)./(wT.^2-w.^2-1i*gT*w);
end

function eps = eps_SiO2_3parOsc(w) % epsilon SiO2 in model of 4-parameter oscillator
    epsinf = 2;
    eps = ones(size(w,1),size(w,2));

    GTOX = [51, 10, 10]; 
    GLOX = [51, 10, 10];
    wTOX = [450, 800, 1045]; 
    wLOX = [505, 830, 1240];   

    for i = 1:length(wTOX)
        eps = eps .* epsilon_1phonon(w, wTOX(i), wLOX(i), GTOX(i), GLOX(i));
    end   
    eps = epsinf * eps;
end