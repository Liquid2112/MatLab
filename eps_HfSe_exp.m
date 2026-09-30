function eps = eps_HfSe_exp(w, Axis) % epsilon 	
	
    switch Axis
            case 'X'                      % 100, x           
                eps = eps_A_GeS_3parOsc(w);  
            case 'Y'                      % 010, y          
                eps = eps_A_GeS_3parOsc(w);  
            case 'Z'                      % 001, z            
                eps = eps_C_GeS_3parOsc(w); 
    end
end

function eps = epsilon_1phonon(w, wT, wL, gT, gL)
    eps  = (wL.^2-w.^2-1i*gL*w)./(wT.^2-w.^2-1i*gT*w);
end

function eps = eps_A_GeS_3parOsc(w) % epsilon  YY GeS in model of 4-parameter oscillator
    epsXinf = 7.25;
    eps = ones(size(w,1),size(w,2));

    GTOX = [6.9]; 
    GLOX = [8.8];
    wTOX = [110.9]; 
    wLOX = [203.8];   

    for i = 1:length(wTOX)
        eps = eps .* epsilon_1phonon(w, wTOX(i), wLOX(i), GTOX(i), GLOX(i));
    end   
    eps = epsXinf * eps;
end

function eps = eps_C_GeS_3parOsc(w) % epsilon ZZ GeS in model of 4-parameter oscillator
    epsZinf = 13.77;
    eps = ones(size(w,1),size(w,2));

    GTOX = [9]; 
    GLOX = [11.9];
    wTOX = [122.1]; 
    wLOX = [141.8];   

    for i = 1:length(wTOX)
        eps = eps .* epsilon_1phonon(w, wTOX(i), wLOX(i), GTOX(i), GLOX(i));
    end 
    eps = epsZinf * eps;
end