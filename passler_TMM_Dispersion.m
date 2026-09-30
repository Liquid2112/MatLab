clear all
close all

Nks = 1; % loop over a parameter, currently thickness
for ks = 1:Nks
    % set up parameters
    d_top    = [30e-9];
    d_bottom = [30e-9];

    % mats   = {'vac','hBN','vac','hBN','CaF2'};     % materials, check epsTarray_generator
    % euler  = {[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]}; % rotation of the layers
    % doping = {[0,0],[0,0],[0,0],[0,0],[0,0]};
    % ds     = {0, d_top(ks), d_bottom(ks), d_top(ks), 0};         % layer thickness

    % mats   = {'vac','hBN','MetaAuC','GaAs'};     % materials, check epsTarray_generator
    % euler  = {[0,0,0],[0,0,0],[0,0,0],[0,0,0]}; % rotation of the layers
    % doping = {[0,0],[0,0],[0,0],[0,0]};
    % ds     = {0, d_top(ks), d_bottom(ks), 0};         % layer thickness

    mats   = {'vac','MetaP5','GaAs'};     % materials, check epsTarray_generator
    euler  = {[0,0,0],[0,0,0],[0,0,0]}; % rotation of the layers
    doping = {[0,0],[0,0],[0,0]};
    ds     = {0, d_top(ks), 0};         % layer thickness

    % mats   = {'vac','MoO3'};     % materials, check epsTarray_generator
    % euler  = {[0,0,0],[0,0,0]}; % rotation of the layers
    % doping = {[0,0],[0,0]};
    % ds     = {0, 0};         % layer thickness

    w  = 400:50:10000;                        % wavenumber, unit: 1/cm
    dz = 0e-9;                              % at which points in z-direction should E be calculated, set 0 if not needed
    % k  = 1.01:0.1:100;                        % relative inplane momentum k>1 means evanescent excitation 

    kn0 = 1.01:100:100000;
    degAngle       = true;   % set true if both in-plane directions are interesting
    doingFFT       = false;   % set false if realspace is not of interest
    analytical     = true;   % set false if the analytical model should not be included in the graph
    ana_harm       = 1;       % how many harmonics should be depicted if analytical is true
    LDOS           = false;

    partial = linspace(0,0,length(w));

    % various declarations
    % Nk   = length(k);
    Nk   = length(kn0);
    Nw   = length(w);
    rppx = zeros(Nk,Nw); rppy = zeros(Nk,Nw);  
    % FFTx = zeros(Nk,Nw); FFTy = zeros(Nk,Nw);  
    
    epsTarray = passler_epsTarray_generator(w,mats,euler,doping);

    l = length(epsTarray(:,1));

    if degAngle == true % permittivity transformation for angle measurement
        if l > 2
            epsH1 = zeros(3,3,Nw,l-2);
            epsH2 = zeros(3,3,Nw,l-2);
            for j=2:l
                epsH1(:,:,:,j-1) = epsTarray{j,1};
                epsH2(:,:,:,j-1) = euler_transform(epsH1(:,:,:,j-1),[0,0,pi/2]);
            end
        elseif l == 2
            epsH1 = epsTarray{2,1};
            epsH2 = euler_transform(epsH1,[0,0,pi/2]);
        end
    end

    for i = 1:Nk            
        % kx0 = k(i).*w;    
        kx0 = kn0(i).*ones(Nw,1); 
        [R,r,t,Ed,zd,qs] = passler_transfer_matrix_modular_k(w,kx0,ds,epsTarray,dz,0,0);
        rppx(i,:) = imag(r(1,:))'; % collect rpp values for momentum kx0 in x-direction
        % FFTx(i,:) = fftshift(ifft(rppx(i,:)));
        if degAngle == true
            if l == 2
                epsTarray{2,1} = epsH2;
                [R,r,t,Ed,zd,qs] = passler_transfer_matrix_modular_k(w,kx0,ds,epsTarray,dz,0,0);
                rppy(i,:) = imag(r(1,:))'; % collect rpp values for momentum kx0 in y-direction  
                % FFTy(i,:) = fftshift(ifft(rppy(i,:)));
                epsTarray{2,1} = epsH1;
            elseif l > 2
                for j=2:l
                    epsTarray{j,1} = epsH2(:,:,:,j-1);
                end
                [R,r,t,Ed,zd,qs] = passler_transfer_matrix_modular_k(w,kx0,ds,epsTarray,dz,0,0);
                rppy(i,:) = imag(r(1,:))'; % collect rpp values for momentum kx0 in y-direction  
                % FFTy(i,:) = fftshift(ifft(rppy(i,:)));
                for j=2:l
                    epsTarray{j,1} = epsH1(:,:,:,j-1);
                end
            end
        end
    end

    if (doingFFT == true) || (LDOS == true)
        if doingFFT == true 
            FFTx = zeros(Nk,Nw); FFTy = zeros(Nk,Nw);  
        end
        if LDOS == true 
            LDOSx = zeros(1,Nw); LDOSy = zeros(1,Nw);  
        end
        for j = 1:Nw
            if doingFFT == true
                FFTx(:,j) = fftshift(ifft(rppx(:,j)));
                if degAngle == true
                    FFTy(:,j) = fftshift(ifft(rppy(:,j)));
                end
            end
            if LDOS == true
                LDOSx(j) = sum(rppx(:,j));
                if degAngle == true
                    LDOSy(j) = sum(rppy(:,j));
                end
            end
        end
    end

    if analytical == true    
        perm{1,1} = epsTarray{1,1};
        perm{3,1} = epsTarray{l,1};

        theta0kreal = zeros(Nw,ana_harm,l-2); 
        theta0kimag = zeros(Nw,ana_harm,l-2); 
        if degAngle == true
            theta90kreal = zeros(Nw,ana_harm,l-2); 
            theta90kimag = zeros(Nw,ana_harm,l-2); 
        end
        for j=2:l-1
            d = ds{1,j};
            perm{2,1} = epsTarray{j,1};
            % [theta0kreal(:,:,j-1), theta0kimag(:,:,j-1)] = dispersionFctnTheory_BNNT(w,d,perm,ana_harm,false,partial);
            [theta0kreal(:,:,j-1), theta0kimag(:,:,j-1)] = dispersionFctnTheory(w,d,perm,ana_harm,false);
            if degAngle == true
                perm{2,1} = epsH2(:,:,:,j-1);
                [theta90kreal(:,:,j-1), theta90kimag(:,:,j-1)] = dispersionFctnTheory(w,d,perm,ana_harm,false);
            end
        end
    end

%% Setting for plots
    
    kh = 0;
    if degAngle == true
        kh = kh + 1;
    end 
    if doingFFT == true
        lhelp = pi/k(1);
        ldash = linspace(-lhelp,lhelp);
        kh = kh + 1;
        if degAngle == true
            kh = kh + 1;
        end
    end 
    if LDOS == true
        kh = kh + 1;
        if degAngle == true
            kh = kh + 1;
        end
    end
    kk = ks*(kh+1) - kh; % running index for figures

%% Plotting

    figure(kk) %
        % imagesc('XData',k,'YData',w,'CData',rppx');
        imagesc('XData',kn0,'YData',w,'CData',rppx');
        xlim([floor(min(k)),floor(max(k))]);
        ylim([floor(min(w)),floor(max(w))]);
        ylabel('Frequency [1/cm]'); xlabel('Momentum kx/k0'); 
        title(['Im(r_{pp})' num2str(ks)]);
        colormap('Jet');
        if (analytical == true) && (l < 6)
            hold on;
            for j=2:l-1
                plot(theta0kreal(:,:,j-1)./(w'*10^2*2*pi),w','--',Color=[1-(j-2)*0.3,1-(j-2)*0.3,0],LineWidth=2);
            end
        end
    
    if doingFFT == true
        kk = kk + 1;
        figure(kk) %
            imagesc('XData',w,'YData',ldash,'CData',abs(FFTx));
            ylim([0,max(ldash)]); 
            xlim([floor(min(w)),floor(max(w))]);
            xlabel('Frequency [1/cm]'); ylabel('Length x [1/cm]'); 
            title(['Im(r_{pp})' num2str(ks)]);
            colormap(parula);
    end

    if LDOS == true
        kk = kk + 1;
        figure(kk)
            plot(w,LDOSx);
    end

    if degAngle == true
        kk = kk + 1;
        figure(kk) %
            % imagesc('XData',k,'YData',w,'CData',rppy');
            imagesc('XData',kn0,'YData',w,'CData',rppy');
            xlim([floor(min(k)),floor(max(k))]);
            ylim([floor(min(w)),floor(max(w))]);
            ylabel('Frequency [1/cm]'); xlabel('Momentum ky/k0'); 
            title(['Im(r_{pp})' num2str(ks)]);
            colormap('Jet');
            if (analytical == true) && (l < 6)
                hold on;
                for j=2:l-1
                    plot(theta90kreal(:,:,j-1)./(w'*10^2*2*pi),w','--',Color=[1-(j-2)*0.3,1-(j-2)*0.3,0],LineWidth=2);
                end
            end
        
        if doingFFT == true
            kk = kk + 1;
            figure(4) %
                imagesc('XData',w,'YData',ldash,'CData',abs(FFTy));
                ylim([0,max(ldash)]);
                xlim([floor(min(w)),floor(max(w))]);
                xlabel('Frequency [1/cm]'); ylabel('Length y [1/cm]'); 
                title(['Im(r_{pp})' num2str(ks)]);
                colormap('Jet');
        end
        if LDOS == true
            kk = kk + 1;
            figure(kk)
                plot(w,LDOSy);
        end
    end
end