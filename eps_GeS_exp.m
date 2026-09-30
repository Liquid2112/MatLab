function eps = eps_GeS_exp(w, Axis) % epsilon 	
	
    switch Axis
            case 'X'                      % 100, x           
                eps = eps_A_GeS_3parOsc(w); 
%                 eps = eps_C_GeS_3parOsc(w);
            case 'Y'                      % 010, y          
                eps = eps_B_GeS_3parOsc(w);  
            case 'Z'                      % 001, z            
                eps = eps_C_GeS_3parOsc(w); 
%                 eps = eps_A_GeS_3parOsc(w); 
    end
end

function eps = epsilon_1phonon(w, wTO, S, gamma)
    eps  = (wTO.^2.*S)./(wTO.^2-w.^2-1i*gamma.*w*wTO);
end

function eps = eps_A_GeS_3parOsc(w) % epsilon XX GeS in model of 4-parameter oscillator

    epsXinf = 13.6;
    eps = zeros(size(w,1),size(w,2));

    gammaX = [.027,.023];
    wTOX = [257.7,119.2]; 
    Sx = [7.6,2.5]; 
%     Sx = [15.2,2.5];

    for i = 1:length(wTOX)
        eps = eps + epsilon_1phonon(w, wTOX(i), Sx(i), gammaX(i));
    end   
    eps = epsXinf + eps; % +40;
end

function eps = eps_B_GeS_3parOsc(w) % epsilon  YY GeS in model of 4-parameter oscillator
    
    epsYinf = 11.8;
    eps = zeros(size(w,1),size(w,2));

    gammaY = [0.022];
    wTOY = [202]; 
    Sy = [17];    % 25

    for i = 1:length(wTOY)
        eps = eps + epsilon_1phonon(w, wTOY(i), Sy(i), gammaY(i));
    end   
    eps = epsYinf + eps;
%     eps = -31*ones(1,length(w));
end

function eps = eps_C_GeS_3parOsc(w) % epsilon ZZ GeS in model of 4-parameter oscillator
    epsZinf = 12;
    eps = zeros(size(w,1),size(w,2));

    gammaZ = [.045,.04 ];    
    wTOZ = [237, 279.7];
    Sz = [7.8, .55];    
    for i = 1:length(wTOZ)
        eps = eps + epsilon_1phonon(w, wTOZ(i), Sz(i), gammaZ(i));
    end  
    eps = epsZinf + eps; %-20
    
%     epsZinf = 12;
%     eps = zeros(size(w,1),size(w,2));
% 
%     gammaZ = [.045];    
%     wTOZ = [237];
%     Sz = [7.8];    
%     for i = 1:length(wTOZ)
%         eps = eps + epsilon_1phonon(w, wTOZ(i), Sz(i), gammaZ(i));
%     end  
%     eps = epsZinf + eps;
end