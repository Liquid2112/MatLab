function [kreal, kimag] = dispersionFctnTheory(wavenumber,thickness,permittivity,ana_harm,values)

    w     = wavenumber;
    perm  = permittivity;
    d     = thickness;   

    c  = 29979245800;
    
    N  = length(w);                        % size of wavenumber vector for calculation
    if values == 1
        q  = zeros(N,ana_harm);
    else
        q  = NaN(N,ana_harm);
    end

    Esup  = perm{1,1};     
    Emat  = perm{2,1};
    Esub  = perm{3,1};
            
    if Emat == Esub
        E1  = Esup(1,1,:);
        Em  = Emat(1,1,:);

        for i = 1:N                         % for whole wavenumber range
            q(i)=sqrt((Em(i)*E1(i)/(Em(i)+E1(i)))*w(i)^2/c^2);
            if (real(q(i))<=0) || (imag(q(i))<=0)
                q(i)=0;          
            end
        end
    else
        E1eff = Esup(1,1,:);
        E1z   = Esup(3,3,:);
        Eeff  = Emat(1,1,:);
        Ez    = Emat(3,3,:);
        E3eff = Esub(1,1,:);
        E3z   = Esub(3,3,:);
        
        for i = 1:N                         % for whole wavenumber range
            rho_l   = 1i*(Ez(i)/Eeff(i))^(1/2);  
            rho_sub = (E3z(i)/E3eff(i))^(1/2); 
            rho_sup = (E1z(i)/E1eff(i))^(1/2); 
            if (real(Eeff(i))>0) && (real(Ez(i))<0)
                Eeff(i) = -Eeff(i);
                Ez(i)   = -Ez(i);
                for j = 1:ana_harm
                    q(i,j) = -rho_l./d*(pi*(j-1)+atan(E1z(i)*rho_l/(Ez(i)*rho_sup))+atan(E3z(i)*rho_l/(Ez(i)*rho_sub)));
                    if (real(q(i,j))<=0) || (real(q(i,j))<=abs(imag(q(i,j))))
                        if values == 1
                            q(i,j) = 0;
                        else
                            q(i,j) = NaN;
                        end
                    end
                end
            elseif (real(Eeff(i))<0) && (real(Ez(i))>0)
                for j = 1:ana_harm
                    q(i,j) = rho_l./d*(pi*(j-1)+atan(E1z(i)*rho_l/(Ez(i)*rho_sup))+atan(E3z(i)*rho_l/(Ez(i)*rho_sub)));
                    if (real(q(i,j))<=0) || (real(q(i,j))<=abs(imag(q(i,j))))
                        if values == 1
                            q(i,j) = 0;
                        else
                            q(i,j) = NaN;
                        end
                    end
                end
            else
                q(i,1) = rho_l./d*(atan(E1z(i)*rho_l/(Ez(i)*rho_sup))+atan(E3z(i)*rho_l/(Ez(i)*rho_sub)));
                if (real(q(i,1))<=0) || (real(q(i,1))<=abs(imag(q(i,1))))
                    if values == 1
                        q(i,1) = 0;
                    else
                        q(i,1) = NaN;
                    end
                end
            end
            % if (real(q(i))<=0) || (imag(q(i))<=0) || (real(q(i))<imag(q(i)))
            
        end
    end
    kreal = real(q);
    kimag = imag(q);
end