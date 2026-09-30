function eps = eps_MoO3_98Mo(w, Axis) % epsilon V2O5 in model of 4-parameter oscillator
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

function eps = epsilon_1phonon(w, wT, wL, gT)
    eps  = (wL.^2-w.^2-1i*gT*w)./(wT.^2-w.^2-1i*gT*w);
end

function eps = eps_A_MoO3_3parOsc(w) % epsilon XX MoO3 in model of 4-parameter oscillator
    epsXinf = 5.3;
    eps = ones(size(w,1),size(w,2));
  
    GTOX = [1.5, 3.8, 49.1]; 
    wTOX = [367, 821, 506.7]; 
    wLOX = [390, 974, 534.3];  

    for i = 1:length(wTOX)
        eps = eps .* epsilon_1phonon(w, wTOX(i), wLOX(i), GTOX(i));
    end   
    eps = epsXinf * eps;
end

function eps = eps_C_MoO3_3parOsc(w) % epsilon  YY MoO3 in model of 4-parameter oscillator
    epsYinf = 10.1;
    eps = ones(size(w,1),size(w,2));
    GTOY = [4, 9.5];
    
    wTOY = [262, 548];
    wLOY = [367, 855];
    for i = 1:length(wTOY)
        eps = eps .* epsilon_1phonon(w, wTOY(i), wLOY(i), GTOY(i));
    end  
    eps = epsYinf * eps;
end

function eps = eps_B_MoO3_3parOsc(w) % epsilon ZZ MoO3 in model of 4-parameter oscillator
    epsZinf = 3.8;
    eps = ones(size(w,1),size(w,2));

    GTOZ = [1, 0.8]; 
    wTOZ = [337, 957];
    wLOZ = [363, 1004];   
    for i = 1:length(wTOZ)
        eps = eps .* epsilon_1phonon(w, wTOZ(i), wLOZ(i), GTOZ(i));
    end  
    eps = epsZinf * eps;
end
