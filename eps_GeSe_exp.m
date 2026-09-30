function eps = eps_GeS_exp(w, Axis) % epsilon 	
	
    switch Axis
            case 'X'                      % 100, x           
                eps = eps_A_GeS_3parOsc(w);  
            case 'Y'                      % 010, y          
                eps = eps_B_GeS_3parOsc(w);  
            case 'Z'                      % 001, z            
                eps = eps_C_GeS_3parOsc(w); 
    end
end

function eps = epsilon_1phonon(w, wTO, S, gamma)
    eps  = (wTO.^2.*S)./(wTO.^2-w.^2-1i*gamma.*w*wTO);
end

function eps = eps_A_GeS_3parOsc(w) % epsilon XX GeSe in model of 4-parameter oscillator

    epsXinf = 13.7;
    eps = zeros(size(w,1),size(w,2));

    gammaX = [.05, .04];
    wTOX = [87.8, 184.5]; 
    Sx = [1.7, 6.62];   

    for i = 1:length(wTOX)
        eps = eps + epsilon_1phonon(w, wTOX(i), Sx(i), gammaX(i));
    end   
    eps = epsXinf + eps;
end

function eps = eps_B_GeS_3parOsc(w) % epsilon  YY GeSe in model of 4-parameter oscillator
    
    epsYinf = 17.5;
    eps = zeros(size(w,1),size(w,2));

    gammaY = [0.054];
    wTOY = [148.7]; 
    Sy = [17.9];    

    for i = 1:length(wTOY)
        eps = eps + epsilon_1phonon(w, wTOY(i), Sy(i), gammaY(i));
    end   
    eps = epsYinf + eps;
end

function eps = eps_C_GeS_3parOsc(w) % epsilon ZZ GeSe in model of 4-parameter oscillator
    epsZinf = 15;
    eps = zeros(size(w,1),size(w,2));

    gammaZ = [.053, .054, 0.03];    
    wTOZ = [83.4, 169.7, 169.9];
    Sz = [1.8, 9.8, 0.43];    
    for i = 1:length(wTOZ)
        eps = eps + epsilon_1phonon(w, wTOZ(i), Sz(i), gammaZ(i));
    end  
    eps = epsZinf + eps;
end