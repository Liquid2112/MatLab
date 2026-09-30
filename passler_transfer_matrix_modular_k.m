function [Rout, rout, tout, Ed, zs, qs] = passler_transfer_matrix_modular_k(w,kx0,ds,epsTarray,dz,dir,fig_num,debug)
        % new version: 15.08.2019, alex 
        % modularize function for better flexibility
        % input arguments: 
        %   w1:        incoming frequency (1/cm)   1d-array
        %   kx0:       inplane momentum   (1/cm)   1d-array
        %   ds:        layer thicknesses  (m)      cell array of doubles, first and last need to be zero
        %   epsTarray: materials          string   cell array of strings, need to match name defined below
        %   dz:        z-resolution for field distributions
        %   dir:       propagation direction  - handle with care. need isotropic, non-absorbing substrate for correct evaluation
        %   fig_num:   figure number for simple debugging, set to zero for no figure output
        %   debug:     high level of debugging, set to zero to avoid a lot of output figures
        % output arguments
        %   Rout:      reflectance
        %              [Rpp, Rss, Rsp, Rps]
        %   rout:      reflection coefficients
        %              [rpp, rps, rss, rsp]
        %   tout:      transmission coefficients/transmitted E-field
        %              amplitudes into the substrate
        %              -> only forward propagating modes in the substrate
        %              [E_fpo, E_fpe, E_fse, E_fso]
        %              each field has x,y,z, components
        %              acronyms:  f: forward,             b: backward
        %                         o: ordinary (or p-pol)  e: extraordinary (or s-pol)
        %                         p: p-pol excitation     s: s-pol excitation
        %   Ed:        cell array (layers) of Efield distributions, same as
        %              above, but additonally backward propgating modes
        %              [E_fpo, E_fpe, E_bpo, E_bpe, E_fso, E_fse, E_bso, E_bse]
        %   zs:        cell array (layers) of z-coordinates for Ed
        %   qs:        cell array (layers) of Berreman q's (kz of the Eigenmodes)
        
    if nargin < 8
        debug = 0;
    end
    % general setting up of the code and declarations

    Cp = @(Ex,Ey)abs(Ex).^2./(abs(Ex).^2+abs(Ey).^2);

    Nds = length(ds);
    Nmats= length(epsTarray);
    Nw = length(w);
    mu1  = 1;

    %% invert the whole structure for backward propagation    
    if dir == -1
        ids = ds;
        iepsTarray = epsTarray;
        for kl = 1:Nds
            ids{kl}        = ds{Nds-kl+1};
            iepsTarray{kl} = epsTarray{Nds-kl+1};
        end
        ds = ids;
        epsTarray = iepsTarray;
    end
    
    
    % set up z-coordinate; 
    % NzS: number of z-points in the substrate;
    NzS = 50;
    z0 = cell(Nds,1); z0{1} = 0;
    z1 = cell(Nds,1); z1{1} = 0;
    zs = cell(Nds,1);
    for kl = 2:Nds
        z0{kl} = z1{kl-1};
        z1{kl} = z0{kl} + ds{kl};
        if kl == Nds
            z1{kl} = z0{kl} + NzS*dz;
        end
        zs{kl} = (z0{kl}):dz:z1{kl};
    end
    z0{1} = - NzS*dz;
    zs{1} = z0{1}:dz:z1{1};


    if Nds ~= Nmats
        fprintf('thickness and material arrays need to be of the same size!\n');
        Rout = []; rout = []; tout = []; Ed = []; zs = []; qs = [];
        return
    end
    M = cell(Nds,1);a = cell(Nds,1);S = cell(Nds,1);
    Dt = cell(Nds,1);qs = cell(Nds,1);  gamma = cell(Nds,1);
    Ai = cell(Nds,1); Ki = cell(Nds,1); Ti = cell(Nds,1); Py = cell(Nds,1);
    psisorted = cell(Nds,1);
    if fig_num
        mymap = zeros(256,3);
        N1 = 100; N2 = 210; 
        for k=1:N1
            mymap(k,:) = [1,(N1-k)/N1,(N1-k)/N1];
        end
        for k=N1:N2
            mymap(k,:) = [1,(k-N1)/(N2-N1),0];
        end
        for k=N2:256
            mymap(k,:) = [1-(k-N2)/(256-N2),1-.5*(k-N2)/(256-N2),0];
        end
    end

    zeros_w   = @(w)zeros(size(w,1),size(w,2));
    ones_w    = @(w)ones(size(w,1),size(w,2));
    epsT_vac = zeros(3,3,Nw);
    epsT_vac(1,:,:) = [ones_w(w);zeros_w(w);zeros_w(w);];      epsT_vac(2,:,:) = [zeros_w(w);ones_w(w);zeros_w(w);];      epsT_vac(3,:,:) = [zeros_w(w);zeros_w(w);ones_w(w);];
    mu1T0 = epsT_vac;

    % determine inplane momentum
    zeta(1,1,:) = kx0(:)./w(:);
    
    % calculate matrices for all layers
    for kl = 1:Nds
        % M matrix
        M{kl} = zeros(6,6,Nw);
        M{kl}(1:3,1:3,:) = epsTarray{kl};
        M{kl}(4:6,4:6,:) = mu1T0;
        % helper in proper dimensions
        d = M{kl}(3,3,:).*M{kl}(6,6,:) - M{kl}(3,6,:).*M{kl}(6,3,:);
        % a matrix    
        a{kl} = zeros(6,6,Nw);
        a{kl}(3,1,:) = (M{kl}(6,1,:).*M{kl}(3,6,:)      - M{kl}(3,1,:).*M{kl}(6,6,:))./d;
        a{kl}(3,2,:) = ((M{kl}(6,2,:)-zeta).*M{kl}(3,6,:)- M{kl}(3,2,:).*M{kl}(6,6,:))./d;
        a{kl}(3,4,:) = (M{kl}(6,4,:).*M{kl}(3,6,:)      - M{kl}(3,4,:).*M{kl}(6,6,:))./d;
        a{kl}(3,5,:) = (M{kl}(6,5,:).*M{kl}(3,6,:)      - (M{kl}(3,5,:)+zeta).*M{kl}(6,6,:))./d;
        a{kl}(6,1,:) = (M{kl}(6,3,:).*M{kl}(3,1,:)      - M{kl}(3,3,:).*M{kl}(6,1,:))./d;
        a{kl}(6,2,:) = (M{kl}(6,3,:).*M{kl}(3,2,:)      - M{kl}(3,3,:).*(M{kl}(6,2,:)-zeta))./d;
        a{kl}(6,4,:) = (M{kl}(6,3,:).*M{kl}(3,4,:)      - M{kl}(3,3,:).*M{kl}(6,4,:))./d;
        a{kl}(6,5,:) = (M{kl}(6,3,:).*(M{kl}(3,5,:)+zeta)- M{kl}(3,3,:).*M{kl}(6,5,:))./d;
        % S matrix
        S{kl} = zeros(4,4,Nw);
        S{kl}(1,1,:) = M{kl}(1,1,:) + M{kl}(1,3,:).*a{kl}(3,1,:) + M{kl}(1,6,:).*a{kl}(6,1,:);
        S{kl}(1,2,:) = M{kl}(1,2,:) + M{kl}(1,3,:).*a{kl}(3,2,:) + M{kl}(1,6,:).*a{kl}(6,2,:);
        S{kl}(1,3,:) = M{kl}(1,4,:) + M{kl}(1,3,:).*a{kl}(3,4,:) + M{kl}(1,6,:).*a{kl}(6,4,:);
        S{kl}(1,4,:) = M{kl}(1,5,:) + M{kl}(1,3,:).*a{kl}(3,5,:) + M{kl}(1,6,:).*a{kl}(6,5,:);
        S{kl}(2,1,:) = M{kl}(2,1,:) + M{kl}(2,3,:).*a{kl}(3,1,:) + (M{kl}(2,6,:)-zeta).*a{kl}(6,1,:);
        S{kl}(2,2,:) = M{kl}(2,2,:) + M{kl}(2,3,:).*a{kl}(3,2,:) + (M{kl}(2,6,:)-zeta).*a{kl}(6,2,:);
        S{kl}(2,3,:) = M{kl}(2,4,:) + M{kl}(2,3,:).*a{kl}(3,4,:) + (M{kl}(2,6,:)-zeta).*a{kl}(6,4,:);
        S{kl}(2,4,:) = M{kl}(2,5,:) + M{kl}(2,3,:).*a{kl}(3,5,:) + (M{kl}(2,6,:)-zeta).*a{kl}(6,5,:);
        S{kl}(3,1,:) = M{kl}(4,1,:) + M{kl}(4,3,:).*a{kl}(3,1,:) + M{kl}(4,6,:).*a{kl}(6,1,:);
        S{kl}(3,2,:) = M{kl}(4,2,:) + M{kl}(4,3,:).*a{kl}(3,2,:) + M{kl}(4,6,:).*a{kl}(6,2,:);
        S{kl}(3,3,:) = M{kl}(4,4,:) + M{kl}(4,3,:).*a{kl}(3,4,:) + M{kl}(4,6,:).*a{kl}(6,4,:);
        S{kl}(3,4,:) = M{kl}(4,5,:) + M{kl}(4,3,:).*a{kl}(3,5,:) + M{kl}(4,6,:).*a{kl}(6,5,:);
        S{kl}(4,1,:) = M{kl}(5,1,:) + (M{kl}(5,3,:)+zeta).*a{kl}(3,1,:) + M{kl}(5,6,:).*a{kl}(6,1,:);
        S{kl}(4,2,:) = M{kl}(5,2,:) + (M{kl}(5,3,:)+zeta).*a{kl}(3,2,:) + M{kl}(5,6,:).*a{kl}(6,2,:);
        S{kl}(4,3,:) = M{kl}(5,4,:) + (M{kl}(5,3,:)+zeta).*a{kl}(3,4,:) + M{kl}(5,6,:).*a{kl}(6,4,:);
        S{kl}(4,4,:) = M{kl}(5,5,:) + (M{kl}(5,3,:)+zeta).*a{kl}(3,5,:) + M{kl}(5,6,:).*a{kl}(6,5,:);
        % Delta matrix
        Dt{kl} = zeros(4,4,Nw);
        Dt{kl}(1,1,:) = S{kl}(4,1,:);  Dt{kl}(1,2,:) = S{kl}(4,4,:);  Dt{kl}(1,3,:) = S{kl}(4,2,:);  Dt{kl}(1,4,:) = - S{kl}(4,3,:);
        Dt{kl}(2,1,:) = S{kl}(1,1,:);  Dt{kl}(2,2,:) = S{kl}(1,4,:);  Dt{kl}(2,3,:) = S{kl}(1,2,:);  Dt{kl}(2,4,:) = - S{kl}(1,3,:);
        Dt{kl}(3,1,:) = -S{kl}(3,1,:); Dt{kl}(3,2,:) = -S{kl}(3,4,:); Dt{kl}(3,3,:) = -S{kl}(3,2,:); Dt{kl}(3,4,:) = S{kl}(3,3,:);
        Dt{kl}(4,1,:) = S{kl}(2,1,:);  Dt{kl}(4,2,:) = S{kl}(2,4,:);  Dt{kl}(4,3,:) = S{kl}(2,2,:);  Dt{kl}(4,4,:) = - S{kl}(2,3,:);

        qs{kl}    = zeros(4,Nw);
        qsd_thr   = 1e-10;
        Py{kl}        = zeros(3,4,Nw);
        psisorted{kl} = zeros(4,4,Nw);
        gamma{kl} = zeros(4,3,Nw);
        Ai{kl}    = zeros(4,4,Nw);
        Ki{kl}    = zeros(4,4,Nw);
        Ti{kl}    = zeros(4,4,Nw);
        locAi = zeros(4,4); locKi = zeros(4,4);

        for kw = 1:Nw
            % determine berremann qi's
            Dtl = zeros(4,4);
            transmode = zeros(1,2);
            reflmode  = zeros(1,2);
            Dtl(:,:) = Dt{kl}(:,:,kw); 
            [psiunsorted, qsunsorted] = eig(Dtl);
            kt = 1; kr = 1;
            % sort berremann qi's
            if any(abs(imag(qsunsorted)))
                for km = 1:4    
                    if imag(qsunsorted(km,km))>=0 
                        transmode(kt) = km; kt = kt + 1;
                    else
                        reflmode(kr) = km; kr = kr +1;
                    end
                end
            else
                for km = 1:4
                    if real(qsunsorted(km,km))>0 
                        transmode(kt) = km; kt = kt + 1;
                    else
                        reflmode(kr) = km; kr = kr +1;
                    end
                end            
            end
            % calculate Poynting vector for each Eigenmode psi
            for km = 1:4
                Ex = psiunsorted(1,km);
                Ey = psiunsorted(3,km);
                Hx = -psiunsorted(4,km);
                Hy = psiunsorted(2,km);
                Ez = a{kl}(3,1,kw)*Ex+a{kl}(3,2,kw)*Ey + a{kl}(3,4,kw)*Hx + a{kl}(3,5,kw)*Hy;
                Hz = a{kl}(6,1,kw)*Ex+a{kl}(6,2,kw)*Ey + a{kl}(6,4,kw)*Hx + a{kl}(6,5,kw)*Hy;
                Pyx = Ey*Hz-Ez*Hy;
                Pyy = Ez*Hx-Ex*Hz;
                Pyz = Ex*Hy-Ey*Hx;
                Py{kl}(1,km,kw)=Pyx;
                Py{kl}(2,km,kw)=Pyy;
                Py{kl}(3,km,kw)=Pyz;
            end
            % check Cp using either the Poynting vector for birefringent
            % materials or the electric field vector for non-birefringent
            % media to sort the modes
            Cp1 = Cp(Py{kl}(1,transmode(1),kw),Py{kl}(2,transmode(1),kw));
            Cp2 = Cp(Py{kl}(1,transmode(2),kw),Py{kl}(2,transmode(2),kw));            
            if abs(Cp1 - Cp2) > qsd_thr % birefringence -> use Poynting vector
                if Cp2 > Cp1
                    transmode([1,2]) = transmode([2,1]);
                end
                if Cp(Py{kl}(1,reflmode(2),kw),Py{kl}(2,reflmode(2),kw)) > Cp(Py{kl}(1,reflmode(1),kw),Py{kl}(2,reflmode(1),kw))
                    reflmode([1,2]) = reflmode([2,1]);
                end
            else                       % no birefringence -> use s-pol/p-pol
                if Cp(psiunsorted(1,transmode(2)),psiunsorted(3,transmode(2))) > Cp(psiunsorted(1,transmode(1)),psiunsorted(3,transmode(1)))
                    transmode([1,2]) = transmode([2,1]);
                end
                if Cp(psiunsorted(1,reflmode(2)),psiunsorted(3,reflmode(2))) > Cp(psiunsorted(1,reflmode(1)),psiunsorted(3,reflmode(1)))
                    reflmode([1,2]) = reflmode([2,1]);
                end
            end
            qs{kl}(:,kw) = [qsunsorted(transmode(1),transmode(1)),qsunsorted(transmode(2),transmode(2)) ...
                            ,qsunsorted(reflmode(1),reflmode(1)),qsunsorted(reflmode(2),reflmode(2))];
            Py{kl}(:,:,kw) = [Py{kl}(:,transmode(1),kw),Py{kl}(:,transmode(2),kw) ...
                            ,Py{kl}(:,reflmode(1),kw),Py{kl}(:,reflmode(2),kw)];
            psisorted{kl}(:,:,kw) = [psiunsorted(:,transmode(1)),psiunsorted(:,transmode(2)) ...
                            ,psiunsorted(:,reflmode(1)),psiunsorted(:,reflmode(2)) ];
            
            % gamma matrix
            gamma{kl}(1,1,:) = ones(Nw,1); gamma{kl}(2,2,:) = ones(Nw,1); 
            gamma{kl}(4,2,:) = ones(Nw,1); gamma{kl}(3,1,:) = -ones(Nw,1);
            % often needed denominator
            den_mueps33mq02 = (mu1*epsTarray{kl}(3,3,kw)-zeta(kw)^2);            
            % check for degenerate qs
            if abs(qs{kl}(1,kw) - qs{kl}(2,kw)) < qsd_thr
                gamma12 = 0.0;
                gamma13 = - (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(1,kw))./den_mueps33mq02;
                gamma21 = 0.0;
                gamma23 = - mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02;
            else            
                gamma12 = (mu1*epsTarray{kl}(2,3,kw)*(mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(1,kw)) ...
                          - mu1*epsTarray{kl}(2,1,kw)*den_mueps33mq02)./(den_mueps33mq02*(mu1*epsTarray{kl}(2,2,kw)-zeta(kw)^2-qs{kl}(1,kw)^2) ...
                          - mu1^2*epsTarray{kl}(2,3,kw)*epsTarray{kl}(3,2,kw));
                if isnan(gamma12), gamma12 = 0.0; end
                gamma13 = - (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(1,kw))./den_mueps33mq02 - mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02*gamma12;
                if isnan(gamma13), gamma13 = - (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(1,kw))./den_mueps33mq02; end
                gamma21 = (mu1*epsTarray{kl}(3,2,kw)*(mu1*epsTarray{kl}(1,3,kw)+zeta(kw)*qs{kl}(2,kw)) - mu1*epsTarray{kl}(1,2,kw)*den_mueps33mq02) ...
                          ./(den_mueps33mq02*(mu1*epsTarray{kl}(1,1,kw)-qs{kl}(2,kw)^2) ...
                          - (mu1*epsTarray{kl}(1,3,kw)+zeta(kw)*qs{kl}(2,kw))*(mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(2,kw)));      
                if isnan(gamma21), gamma21 = 0.0; end
                gamma23 = - (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(2,kw))./den_mueps33mq02*gamma21 - mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02;
                if isnan(gamma23), gamma23 = - mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02; end
            end    
            if abs(qs{kl}(3,kw) - qs{kl}(4,kw)) < qsd_thr
                gamma32 = 0.0;
                gamma33 = (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(3,kw))./den_mueps33mq02;
                gamma41 = 0.0;
                gamma43 = - mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02;
            else            
                gamma32 = (mu1*epsTarray{kl}(2,1,kw)*den_mueps33mq02 - mu1*epsTarray{kl}(2,3,kw)*(mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(3,kw))) ...
                          ./(den_mueps33mq02*(mu1*epsTarray{kl}(2,2,kw)-zeta(kw)^2-qs{kl}(3,kw)^2) - mu1^2*epsTarray{kl}(2,3,kw)*epsTarray{kl}(3,2,kw));
                if isnan(gamma32), gamma32 = 0.0; end
                gamma33 = (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(3,kw))./den_mueps33mq02 + mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02*gamma32;
                if isnan(gamma33), gamma33 = (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(3,kw))./den_mueps33mq02; end
                gamma41 = (mu1*epsTarray{kl}(3,2,kw)*(mu1*epsTarray{kl}(1,3,kw)+zeta(kw)*qs{kl}(4,kw)) - mu1*epsTarray{kl}(1,2,kw)*den_mueps33mq02) ...
                          ./(den_mueps33mq02*(mu1*epsTarray{kl}(1,1,kw)-qs{kl}(4,kw)^2) ...
                          - (mu1*epsTarray{kl}(1,3,kw)+zeta(kw)*qs{kl}(4,kw))*(mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(4,kw)));     
                if isnan(gamma41), gamma41 = 0.0; end
                gamma43 = - (mu1*epsTarray{kl}(3,1,kw)+zeta(kw)*qs{kl}(4,kw))./den_mueps33mq02*gamma41 - mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02;
                if isnan(gamma43), gamma43 =  - mu1*epsTarray{kl}(3,2,kw)./den_mueps33mq02; end
            end
            
            gamma{kl}(1,2,kw) = gamma12;  gamma{kl}(1,3,kw) = gamma13;   
            gamma{kl}(2,1,kw) = gamma21;  gamma{kl}(2,3,kw) = gamma23;   
            gamma{kl}(3,2,kw) = gamma32;  gamma{kl}(3,3,kw) = gamma33;
            gamma{kl}(4,1,kw) = gamma41;  gamma{kl}(4,3,kw) = gamma43;

            % normalize gamma vectors
            % Max solution
            for km = 1:4
                if max(isinf(real(gamma{kl}(km,:,kw))))==1
                    gamma{kl}(km,:,kw)=1;
                    % disp(1);
                else
                    gamma{kl}(km,:,kw)=gamma{kl}(km,:,kw)/norm(gamma{kl}(km,:,kw));
                    % disp(2);
                end
            end
            
            % transfer matrices for each layer
            Ai{kl}(1:2,1:4,kw) = permute(gamma{kl}(1:4,1:2,kw),[2,1,3]);
            for k = 1:4 
                Ai{kl}(3,k,kw) = 1/mu1*(qs{kl}(k,kw).*gamma{kl}(k,1,kw) - zeta(kw)*gamma{kl}(k,3,kw)); 
                Ai{kl}(4,k,kw) = 1/mu1*qs{kl}(k,kw).*gamma{kl}(k,2,kw);
                Ki{kl}(k,k,kw) = exp(-1i*2*pi*w(kw)*1e2*qs{kl}(k,kw)*ds{kl});
            end

            locAi(:,:) = Ai{kl}(:,:,kw); locKi(:,:) = Ki{kl}(:,:,kw); 
            Ti{kl}(:,:,kw) = locAi*locKi/locAi;        
        end
    end
    
    % calculate total transfer matrix
    A0 = Ai{1};
    Af = Ai{Nds};
    SwapColumns1324=[[1,0,0,0];[0,0,1,0];[0,1,0,0];[0,0,0,1];];
    unity44 = [[1,0,0,0];[0,1,0,0];[0,0,1,0];[0,0,0,1];];
    T = zeros(4,4,Nw);
    Tl = zeros(4,4); Tsl = zeros(4,4); A0l = zeros(4,4); Afl = zeros(4,4);

    for k = 1:4
        T(k,k,:) = 1.0;
    end
    for kw = 1:Nw
        Tl(:,:)  = T(:,:,kw);
        for kl = (Nds-1):-1:2
            Tsl(:,:) = Ti{kl}(:,:,kw);
            Tl = Tsl*Tl;
        end
        A0l(:,:) = A0(:,:,kw); Afl(:,:) = Af(:,:,kw);
        Tl = unity44/A0l*Tl*Afl;
        T(:,:,kw) = SwapColumns1324*Tl*SwapColumns1324;
    end
    % disp(T);

    %% evaluate reflectance, reflection and transmission coefficients
    rM = zeros(2,2,Nw);
    rMden = 1./(T(1,1,:).*T(3,3,:) - T(1,3,:).*T(3,1,:));
    rM(1,1,:) = (T(2,1,:).*T(3,3,:) - T(2,3,:).*T(3,1,:)).*rMden;
    rM(1,2,:) = (T(4,1,:).*T(3,3,:) - T(4,3,:).*T(3,1,:)).*rMden;
    rM(2,1,:) = (T(1,1,:).*T(2,3,:) - T(2,1,:).*T(1,3,:)).*rMden;
    rM(2,2,:) = (T(1,1,:).*T(4,3,:) - T(4,1,:).*T(1,3,:)).*rMden;

    tM = zeros(2,2,Nw);
    tM(1,1,:) = T(3,3,:).*rMden;
    tM(1,2,:) = -T(3,1,:).*rMden;
    tM(2,1,:) = -T(1,3,:).*rMden;
    tM(2,2,:) = T(1,1,:).*rMden;
    
    % reflectance
    Rout = zeros(6,Nw);
    Rout(1,:) = abs(permute(rM(1,1,:).^2,[1,3,2])); % ppol reflectance 
    Rout(2,:) = abs(permute(rM(2,2,:).^2,[1,3,2])); % spol reflectance
    Rout(3,:) = abs(permute(rM(2,1,:).^2,[1,3,2])); % s->p pol reflectance
    Rout(4,:) = abs(permute(rM(1,2,:).^2,[1,3,2])); % p->s spol reflectance
    
    % reflection coefficients (field)
    % since the incident medium is isotropic and non-absorbing, it is
    % sufficient to evaluate the mode amplitudes (which are equivalent to
    % the x-field (ppol) and y-field (spol) amplitudes here.
    rout(1,:) = permute(rM(1,1,:),[1,3,2]); % ppol reflection coefficient
    rout(2,:) = permute(rM(1,2,:),[1,3,2]); % p->spol reflection coefficient
    rout(3,:) = permute(rM(2,2,:),[1,3,2]); % spol reflection coefficient
    rout(4,:) = permute(rM(2,1,:),[1,3,2]); % s->ppol reflection coefficient    

    % transmission coefficients (field)
    % since the substrate can be birefringent, we output the full
    % x-,y-,z-fields for both modes (e/o,s/p) seperately 
    % to calculate Lzz, one needs to additionally account for the incidence
    % angle (ratio of Ex and Ez in the incidence medium)
    toutpp = zeros(3,Nw);toutps = zeros(3,Nw);toutss = zeros(3,Nw);toutsp = zeros(3,Nw);
    % p-pol incidence
    % ppol->ordinary/ppol->ppol transmission coefficients
    toutpp(1,:) = permute(tM(1,1,:).*gamma{end}(1,1,:),[1,3,2]);
    toutpp(2,:) = permute(tM(1,1,:).*gamma{end}(1,2,:),[1,3,2]);
    toutpp(3,:) = permute(tM(1,1,:).*gamma{end}(1,3,:),[1,3,2]);        
    % ppol->extraordinary/ppol->spol transmission coefficients
    toutps(1,:) = permute(tM(1,2,:).*gamma{end}(2,1,:),[1,3,2]);
    toutps(2,:) = permute(tM(1,2,:).*gamma{end}(2,2,:),[1,3,2]);
    toutps(3,:) = permute(tM(1,2,:).*gamma{end}(2,3,:),[1,3,2]);        
    % s-pol incidence
    % spol->extraordinary/spol->spol transmission coefficients
    toutss(1,:) = permute(tM(2,2,:).*gamma{end}(2,1,:),[1,3,2]);
    toutss(2,:) = permute(tM(2,2,:).*gamma{end}(2,2,:),[1,3,2]);
    toutss(3,:) = permute(tM(2,2,:).*gamma{end}(2,3,:),[1,3,2]);
    % spol->ordinary/spol->ppol transmission coefficients
    toutsp(1,:) = permute(tM(2,1,:).*gamma{end}(1,1,:),[1,3,2]);
    toutsp(2,:) = permute(tM(2,1,:).*gamma{end}(1,2,:),[1,3,2]);
    toutsp(3,:) = permute(tM(2,1,:).*gamma{end}(1,3,:),[1,3,2]);
    
    tout = [toutpp;toutps;toutss;toutsp;];
    
    if fig_num
        figure(fig_num), 
        subplot(2,2,1), hold off, plot(w,real(permute(rM(1,1,:),[3,1,2])))
        hold on, plot(w,imag(permute(rM(2,1,:),[3,1,2]))) 
        subplot(2,2,3), plot(w,Rout(1,:),w,Rout(3,:)), ylim([0,1])
        subplot(2,2,2), hold off, 
        plot(w,real(tout),w,imag(tout)),
        subplot(2,2,4), hold off, 
        plot(w,abs(tout)),set(gca,'YScale','log');
        xlim([min(w),max(w)]); 
    end
        
    %% calc field distributions
    % first: calculate and propagation basis vector amplitudes for p-pol
    % and s-pol excitation
    % second: project onto the basis vectors in x-,y-,z- basis
    Ed = cell(Nds,1);                    % field distributions in each layer
    E0  = cell(Nds,1); E1 = cell(Nds,1); % fields at the front and back interface for each layer    
    for kl = Nds:-1:1
        Nz = length(zs{kl});
        Ed{kl} = zeros(24,Nw,Nz);
        E0{kl} = zeros(8,Nw); E1{kl} = zeros(8,Nw);
        locKi= zeros(4,4); 
        % first calculate all interface field amplitudes using the transfer matrices 
        if kl == Nds
            % first four components: p-pol excitation
            E0{kl}(1,:) = permute(tM(1,1,:),[2,3,1]);
            E0{kl}(2,:) = permute(tM(1,2,:),[2,3,1]);
            % second four components: s-pol excitation
            E0{kl}(5,:) = permute(tM(2,1,:),[2,3,1]);
            E0{kl}(6,:) = permute(tM(2,2,:),[2,3,1]);            
        else
            for kw = 1:Nw
                Aim1l = Ai{kl}(:,:,kw);
                Ail   = Ai{kl+1}(:,:,kw);
                Li   = unity44/Aim1l*Ail;
                E0l(:,1)  = E0{kl+1}(:,kw); % Efield at the interface in layer N+1;
                % Efield on the other side of that interface, seperately
                % for p-pol excitation and s-pol excitation
                E1l(1:4,1)  = Li*E0l(1:4,1); 
                E1l(5:8,1)  = Li*E0l(5:8,1); 
                E1{kl}(:,kw) = E1l;
                locKi(:,:) = Ki{kl}(:,:,kw); % propgation back through layer N
                if kl == 1
                    for k = 1:4
                        locKi(k,k) = exp(-1i*2*pi*w(kw)*1e2*qs{kl}(k,kw)*(z1{kl}-z0{kl}));
                    end
                end
                % Efield at the front of layer N, 
                % for s- and p-pol excitation
                E0{kl}(1:4,kw) = locKi*E1l(1:4); 
                E0{kl}(5:8,kw) = locKi*E1l(5:8);       
            end
        end
        dKiz = zeros(4,4);
        % now propagate, always front surface (E0) to back surface (E1)
        % in each layer, conjugate propagator from the main matrix
        for kw = 1:Nw
            for k = 1:4
                dKiz(k,k) = exp(-1i*2*pi*w(kw)*1e2*qs{kl}(k,kw)*(-dz));
            end
            E0l(:,1) = E0{kl}(:,kw);
            Edl = E0l;
            for kz = 1:Nz
                % collect all fields components
                % fields for p-pol incidence
                    % forward, ordinary/ppol 
                Ed{kl}(1,kw,kz) = Edl(1,1).*gamma{kl}(1,1,kw);
                Ed{kl}(2,kw,kz) = Edl(1,1).*gamma{kl}(1,2,kw);
                Ed{kl}(3,kw,kz) = Edl(1,1).*gamma{kl}(1,3,kw);
                    % forward, extraordinary/spol 
                Ed{kl}(4,kw,kz) = Edl(2,1).*gamma{kl}(2,1,kw);
                Ed{kl}(5,kw,kz) = Edl(2,1).*gamma{kl}(2,2,kw);
                Ed{kl}(6,kw,kz) = Edl(2,1).*gamma{kl}(2,3,kw);                
                    % backward, ordinary/ppol
                Ed{kl}(7,kw,kz) = Edl(3,1).*gamma{kl}(3,1,kw);
                Ed{kl}(8,kw,kz) = Edl(3,1).*gamma{kl}(3,2,kw);
                Ed{kl}(9,kw,kz) = Edl(3,1).*gamma{kl}(3,3,kw);
                    % backward, extraordinary/spol
                Ed{kl}(10,kw,kz) = Edl(4,1).*gamma{kl}(4,1,kw);
                Ed{kl}(11,kw,kz) = Edl(4,1).*gamma{kl}(4,2,kw);
                Ed{kl}(12,kw,kz) = Edl(4,1).*gamma{kl}(4,3,kw);                
                % fields for s-pol incidence
                    % forward, ordinary/ppol
                Ed{kl}(13,kw,kz) = Edl(5,1).*gamma{kl}(1,1,kw);
                Ed{kl}(14,kw,kz) = Edl(5,1).*gamma{kl}(1,2,kw);  
                Ed{kl}(15,kw,kz) = Edl(5,1).*gamma{kl}(1,3,kw);
                    % forward, extraordinary/spol
                Ed{kl}(16,kw,kz) = Edl(6,1).*gamma{kl}(2,1,kw);
                Ed{kl}(17,kw,kz) = Edl(6,1).*gamma{kl}(2,2,kw);
                Ed{kl}(18,kw,kz) = Edl(6,1).*gamma{kl}(2,3,kw);                
                    % backward, ordinary/ppol
                Ed{kl}(19,kw,kz) = Edl(7,1).*gamma{kl}(3,1,kw);
                Ed{kl}(20,kw,kz) = Edl(7,1).*gamma{kl}(3,2,kw);
                Ed{kl}(21,kw,kz) = Edl(7,1).*gamma{kl}(3,3,kw);
                    % backward, extraordinary/spol
                Ed{kl}(22,kw,kz) = Edl(8,1).*gamma{kl}(4,1,kw);
                Ed{kl}(23,kw,kz) = Edl(8,1).*gamma{kl}(4,2,kw);
                Ed{kl}(24,kw,kz) = Edl(8,1).*gamma{kl}(4,3,kw);                                    
                
                % now propagate mode amplitudes at z and w in layer N                    
                Edl(1:4,1) = dKiz*Edl(1:4,1);   % p-pol excitation
                Edl(5:8,1) = dKiz*Edl(5:8,1);   % s-pol excitation                
            end
        end                
        if fig_num && kl>1
            figure(667), subplot(Nds-1,1,kl-1), imagesc(w,zs{kl},permute(abs(Ed{kl}(1,:,:)-Ed{kl}(4,:,:)),[3,2,1]),[0,8])
        end
    end
    %now invert the coordinate system for backward propagation
    if dir == -1
        iEd = Ed;
        izs = zs;
        iqs = qs;
        for kl = 1:Nds
            iEd{kl} =  Ed{Nds-kl+1};
            izs{kl} = -zs{Nds-kl+1};
            iqs{kl} =  qs{Nds-kl+1};
        end
        Ed = iEd;
        zs = izs;
        qs = iqs;
    end
end