function eps = eps_MoO3_Full(w, Axis) % epsilon V2O5 in model of 4-parameter oscillator
    % data from Gonzalo
    switch Axis
            case 'X'                      % 100              
                eps = eps_A_MoO3_3parOsc(w);  
            case 'Z'                      % 010             
                eps = eps_B_MoO3_3parOsc(w);  
            case 'Y'                      % 001             
                eps = eps_C_MoO3_3parOsc(w); 
    end
end

function eps = epsilon_1phonon(w, wT, wL, gT, gL)
    eps  = (wL.^2-w.^2-1i*gL*w)./(wT.^2-w.^2-1i*gT*w);
end

function eps = eps_A_MoO3_3parOsc(w) % epsilon XX MoO3 in model of 4-parameter oscillator
    epsXinf = 5.78;
    eps = ones(size(w,1),size(w,2));

    GTOX = [3, 6, .35, 49.1]; 
    GLOX = [3, 6, .35, 49.1];
    wTOX = [367, 821.4, 998.7, 506.7]; 
    wLOX = [390, 963.0, 999.2, 534.3];   
    % GTOX = [1.5, 6, 49.1]; 
    % GLOX = [1.5, 6, 49.1];
    % wTOX = [367, 821.4, 506.7]; 
    % wLOX = [390, 963.0, 534.3];  

    for i = 1:length(wTOX)
        eps = eps .* epsilon_1phonon(w, wTOX(i), wLOX(i), GTOX(i), GLOX(i));
    end   
    eps = epsXinf * eps;
end

function eps = eps_C_MoO3_3parOsc(w) % epsilon  YY MoO3 in model of 4-parameter oscillator
    epsYinf = 6.07;
    eps = ones(size(w,1),size(w,2));
    GTOY = [4, 9.5];
    GLOY = [4, 9.5];
    
    wTOY = [262, 544.6];
    wLOY = [367, 850.1];
    for i = 1:length(wTOY)
        eps = eps .* epsilon_1phonon(w, wTOY(i), wLOY(i), GTOY(i), GLOY(i));
    end  
    eps = epsYinf * eps;
end

function eps = eps_B_MoO3_3parOsc(w) % epsilon ZZ MoO3 in model of 4-parameter oscillator
    epsZinf = 4.47;
    eps = ones(size(w,1),size(w,2));

    GTOZ = [1, 1.5];
    GLOZ = [1, 1.5];  
    wTOZ = [337, 956.7];
    wLOZ = [363, 1006.9];   
%     wTOZ = [307, 956.7];
%     wLOZ = [333, 1006.9]; 
    for i = 1:length(wTOZ)
        eps = eps .* epsilon_1phonon(w, wTOZ(i), wLOZ(i), GTOZ(i), GLOZ(i));
    end  
    eps = epsZinf * eps;
end