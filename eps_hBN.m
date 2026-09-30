function eps = eps_hBN(w, Axis) % epsilon V2O5 in model of 4-parameter oscillator
    % data from Gonzalo
    switch Axis
            case 'ip'                      % 100              
                eps = eps_ip_BNNT_3parOsc(w);  
            case 'oop'                      % 010             
                eps = eps_oop_BNNT_3parOsc(w);  
    end
end

% function eps = epsilon_1phonon(w, wT, wL, gamma)
%     eps  = (wL.^2-wT.^2)./(wT.^2-w.^2-1i*gT*w);
% end

function eps = eps_ip_BNNT_3parOsc(w) % epsilon XX MoO3 in model of 4-parameter oscillator  
 
    gamma = [7]; % matrix generator
    wTOip  = [1360]; 
    wLOip  = [1614];  
    epsIPinf = 4.9;

    % p = partial;
    
    % gamma = [4]; % BNNT publication
    % wTOip  = [1362]; 
    % wLOip  = [1527]; 
    % epsIPinf = 4.95;

    % for i = 1:length(wTOX)
    %     eps = eps .* epsilon_1phonon(w, wTOX(i), wLOX(i), gamma(i));
    % end   
    % eps = epsIPinf * eps;
    eps = epsIPinf .* (1 + (wLOip.^2 - wTOip.^2)./(wTOip.^2 - w.^2 - 1i*w*gamma));
    % if partial ~= 0
    %     % p = zeros(length(w),1)';
    %     % for i=1:length(w)
    %     %     if (w(i)==wLOip)
    %     %         LOp = i;
        %     end
    %     %     if (w(i)==1450) %(w(i)==wTOip)
    %     %         TOp = i;
    %     %     end
    %     % end
    %     % lRB = LOp-TOp;
    %     % p(TOp:LOp) = linspace(partial,0,lRB+1)';
    %     eps = (eps.*(2.*p.*(1-eps)+1+2*eps))./(2*eps+1-p.*(1-eps));
    % end
end

function eps = eps_oop_BNNT_3parOsc(w) % epsilon  YY MoO3 in model of 4-parameter oscillator
    epsOOPinf = 4.1;
    % eps = ones(size(w,1),size(w,2));

    gamma = [1];
    wTOoop  = [760]; 
    wLOoop  = [811];  

    % for i = 1:length(wTOY)
    %     eps = eps .* epsilon_1phonon(w, wTOY(i), wLOY(i), gamma(i));
    % end  
    % eps = epsOOPinf * eps;
    eps = epsOOPinf .* (1 + (wLOoop.^2 - wTOoop.^2)./(wTOoop.^2 - w.^2 - 1i*w*gamma));
    % eps = (eps.*(2*partial.*(1-eps)+1+2*eps))./(2*eps+1-partial.*(1-eps));
end