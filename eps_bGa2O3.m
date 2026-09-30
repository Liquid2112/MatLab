function epsT = eps_bGa2O3(w,N)

debug = 0;

Nw = length(w);

if N(1)=='T'
    if N(2) == 'i'
        %data = importdata(['e',N(3:6),' - 296 isotropic mobility.txt']);
        data = importdata(['bGOepsN',N(3:6),'i.txt']);
        data = data.data;
    elseif N(2) == 'a'
        data = importdata(['bGOepsN',N(3:6),'a296-37-37.txt']);
        data = data.data;
    elseif N(2) == 'r'
        data = importdata(['e',N(3:6),' - 37-37-296 mobility.txt']);
        data = data.data;
    elseif N(2) == 'm'
        data = importdata(['e',N(3:6),' - 296 isotropic mobility.txt']);
        data = data.data;        
    elseif N(2) == 'u'
        data = importdata('Undoped (theoretical).txt');
    end
else
    switch N
        case '0'
            load('epsparms_bGa2O3.mat');
        case '1' 
            load('epsparms_bGa2O3_1.mat');
        case '2' 
            load('epsparms_bGa2O3_2.mat');
        case '3' 
            load('epsparms_bGa2O3_3.mat');
        case '12' 
            load('epsparms_bGa2O3_12.mat');
        case '13' 
            load('epsparms_bGa2O3_13.mat');
        case '23' 
            load('epsparms_bGa2O3_23.mat');
        case '123' 
            load('epsparms_bGa2O3_123.mat');
        case '123i'
            load('epsparms_bGa2O3_123o.mat');    
        case 'f0'
            load('epsparms_bGa2O3_f0.mat');
        case 'f0.2'
            load('epsparms_bGa2O3_f0.2.mat');
        case 'f0.4'
            load('epsparms_bGa2O3_f0.4.mat');
        case 'f0.6'
            load('epsparms_bGa2O3_f0.6.mat');            
        case 'f0.8'
            load('epsparms_bGa2O3_f0.8.mat');
        case 'f1'
            load('epsparms_bGa2O3_f1.mat');
        case 'VWASE'
            load('epsparms_bGO_VWASE.mat');
        case 'lossless'
            load('epsparms_bGa2O3_lossless.mat');
        case 'N0'
            data = importdata('bGOepsN0.txt');
        case 'N17'
            data = importdata('bGOepsN17.txt');
        case 'N18'
            data = importdata('bGOepsN18.txt');
        case 'N19'
            data = importdata('bGOepsN19.txt');
        case 'N18.0i'
            data = importdata('bGOepsN18.0i.txt');
        case 'N18.1i'
            data = importdata('bGOepsN18.1i.txt');
        case 'N18.2i'
            data = importdata('bGOepsN18.2i.txt');
        case 'N18.3i'
            data = importdata('bGOepsN18.3i.txt');
        case 'N18.4i'
            data = importdata('bGOepsN18.4i.txt');
        case 'N18.5i'
            data = importdata('bGOepsN18.5i.txt');
        case 'N18.6i'
            data = importdata('bGOepsN18.6i.txt');
        case 'N18.6a-148-37-37'
            data = importdata('bGOepsN18.6a148-37-37.txt');
        case 'N18.6a-37-37-148'
            data = importdata('bGOepsN18.6a37-37-148.txt'); 
        case 'N18.6a-74-37-37'
            data = importdata('bGOepsN18.6a74-37-37.txt');
        case 'N18.6a-37-37-74'
            data = importdata('bGOepsN18.6a37-37-74.txt');                  
    end
end
if (N(1) == 'N') || (N(1) == 'T') %tabulated data
    if N(1) == 'T'
        epsdata = data(2:end,2:end);
        wdata   = data(2:end,1);
    else    
        epsdata = data.data(2:end,2:end);
        wdata   = data.data(2:end,1);
    end
    eps_loc = cell(9,1);
    for kk = 1:2:17
        eps_loc{(kk+1)/2} = interp1(wdata',epsdata(:,kk)'+1i*epsdata(:,kk+1)',w);
    end
    epsT = zeros(3,3,Nw);
    for kk = 1:3
        for ii = 1:3
            epsT(kk,ii,:) = eps_loc{3*(kk-1)+ii};
        end
    end
elseif or(N(1) == 'V',N(1) == 'f')
    epsinf = epsinf_bGa2O3;
    epsparms = epsparms_bGa2O3.Variables;
    
    epsparms(:,6) = pi/180*epsparms(:,6);
    epsparms(:,7) = pi/180*epsparms(:,7);
    epsT = zeros(3,3,Nw);
    for kk = 1:3
        for jj = 1:3
            epsT(kk,jj,:) = epsinf(kk,jj);
        end
    end

    for ii = 1:length(epsparms(:,1))
        osc = permute((epsparms(ii,3).^2+1i*epsparms(ii,5).*w)./(epsparms(ii,2).^2 - w.^2 - 1i*epsparms(ii,4).*w),[1,3,2]);
        epsT(1,1,:) = epsT(1,1,:) + sin(epsparms(ii,7)).^2.*cos(epsparms(ii,6)).^2.*osc;
        epsT(2,2,:) = epsT(2,2,:) + sin(epsparms(ii,7)).^2.*sin(epsparms(ii,6)).^2.*osc;
        epsT(3,3,:) = epsT(3,3,:) + cos(epsparms(ii,7)).^2.*osc;
        epsT(1,2,:) = epsT(1,2,:) + sin(epsparms(ii,7)).^2.*sin(epsparms(ii,6)).*cos(epsparms(ii,6)).*osc;
        epsT(1,3,:) = epsT(1,3,:) + sin(epsparms(ii,7)).*cos(epsparms(ii,7)).*cos(epsparms(ii,6)).*osc;
        epsT(2,3,:) = epsT(2,3,:) + sin(epsparms(ii,7)).*cos(epsparms(ii,7)).*sin(epsparms(ii,6)).*osc;
        if debug
            fprintf('\n ii = %d',ii);
            epsT(:,:,4)
            figure(1), 
            %subplot(2,1,1), plot(w,permute(real(osc),[3,1,2]),w,permute(imag(osc),[3,1,2]))
            for kk = 1:3
            subplot(2,3,kk), plot(w,permute(real(epsT(kk,kk,:)),[3,1,2]))
            subplot(2,3,3+kk), plot(w,permute(imag(epsT(kk,kk,:)),[3,1,2]))
            %set(gca,'XDir','reverse')
            end
        end
    end

    epsT(2,1,:) = epsT(1,2,:);
    epsT(3,1,:) = epsT(1,3,:);
    epsT(3,2,:) = epsT(2,3,:);  
else% model data
    epsinf = epsinf_bGa2O3;
    epsparms = epsparms_bGa2O3.Variables;

    epsparms(:,5) = pi/180*epsparms(:,5);
    epsparms(:,6) = pi/180*epsparms(:,6);
    
    epsT = zeros(3,3,Nw);
    for kk = 1:3
        for jj = 1:3
            epsT(kk,jj,:) = epsinf(kk,jj);
        end
    end

    for ii = 1:length(epsparms(:,1))
        osc = permute(epsparms(ii,3).^2./(epsparms(ii,2).^2 - w.^2 - 1i*epsparms(ii,4).*w),[1,3,2]);
        epsT(1,1,:) = epsT(1,1,:) + sin(epsparms(ii,6)).^2.*cos(epsparms(ii,5)).^2.*osc;
        epsT(2,2,:) = epsT(2,2,:) + sin(epsparms(ii,6)).^2.*sin(epsparms(ii,5)).^2.*osc;
        epsT(3,3,:) = epsT(3,3,:) + cos(epsparms(ii,6)).^2.*osc;
        epsT(1,2,:) = epsT(1,2,:) + sin(epsparms(ii,6)).^2.*sin(epsparms(ii,5)).*cos(epsparms(ii,5)).*osc;
        epsT(1,3,:) = epsT(1,3,:) + sin(epsparms(ii,6)).*cos(epsparms(ii,6)).*cos(epsparms(ii,5)).*osc;
        epsT(2,3,:) = epsT(2,3,:) + sin(epsparms(ii,6)).*cos(epsparms(ii,6)).*sin(epsparms(ii,5)).*osc;
        if debug
            fprintf('\n ii = %d',ii);
            epsT(:,:,4)
            figure(1), 
            %subplot(2,1,1), plot(w,permute(real(osc),[3,1,2]),w,permute(imag(osc),[3,1,2]))
            for kk = 1:3
            subplot(2,3,kk), plot(w,permute(real(epsT(kk,kk,:)),[3,1,2]))
            subplot(2,3,3+kk), plot(w,permute(imag(epsT(kk,kk,:)),[3,1,2]))
            %set(gca,'XDir','reverse')
            end
        end
    end

    epsT(2,1,:) = epsT(1,2,:);
    epsT(3,1,:) = epsT(1,3,:);
    epsT(3,2,:) = epsT(2,3,:);
end

if debug
    f2 = figure(2);
    close(f2);
    figure(2), 
    %subplot(2,1,1), plot(w,permute(real(osc),[3,1,2]),w,permute(imag(osc),[3,1,2]))
    for kk = 1:3
        subplot(2,2,1), hold on, plot(w,permute(real(epsT(kk,kk,:)),[3,1,2]))
        xlim([850,980]);
        set(gca,'XDir','reverse')
        subplot(2,2,2), hold on, plot(w,permute(real(epsT(kk,kk,:)),[3,1,2]))
        xlim([0,850]);        
        set(gca,'XDir','reverse')        
        subplot(2,2,3), hold on, plot(w,permute(imag(epsT(kk,kk,:)),[3,1,2]))
        xlim([850,980]);        
        set(gca,'XDir','reverse')        
        subplot(2,2,4), hold on, plot(w,permute(imag(epsT(kk,kk,:)),[3,1,2]))
        xlim([0,850]);
        set(gca,'XDir','reverse')                
    end
    legend('\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}');
end


end

