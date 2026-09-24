function controls = getControls(type_of_control)
    % type_of_control can be 'pi' 'alphaRIO' and 'alphadeltaRIO'
    if strcmpi(type_of_control,'alphadeltaRIO')
        global  lam1 lam2 lam3 lam4
        lam1 =  0.208123297235226; lam2 =  -0.388787906696785;
        lam3 =  0.320148587892857; lam4 =  -0.598058460777178;
        gami=pi/2;  thi=0;
        dgami=0;    dthi=1;
        T_max = 6.69366153065;
        t = linspace(0, T_max, 20001);
        
        options=odeset('RelTol',1e-10,'AbsTol',1e-10);
        [t,GT]=ode45('odefun_alphadeltaRIO',t,[gami;thi+1e-14;dgami;dthi],options);
        gam=GT(:,1);   th=GT(:,2);  
        dgam=GT(:,3);  dth=GT(:,4);
        Omega=sqrt(dth.^2+dgam.^2.*sin(th).^2);
        Delta=sin(th).*(lam1*sin(gam)-lam2*cos(gam))+(dth.*(lam3*sin(gam)-lam4*cos(gam))+dgam.*cos(th).*sin(th).*(lam3*cos(gam)+lam4*sin(gam)))./Omega.^2;
        controls =[t,Delta,Omega];

    elseif strcmpi(type_of_control,'alphaRIO') 
        Area = 5.841909888;
        Omega0 = 1;
        T_max = Area/Omega0;
        t = linspace(0, T_max, 20001);
    
        Omega=0.*t+Omega0; % constant pulse
        m=0.234883370; 
        D0=8*ellipticK(m)*sqrt(m)/T_max;
        Delta=D0*jacobiCN(4*ellipticK(m).*t/T_max+ellipticK(m),m);
        controls =[t',Delta',Omega'];
    
    elseif strcmpi(type_of_control,'pi')
        Area = pi;
        Omega0=1; 
        T_max = Area/Omega0;
        t = linspace(0, T_max, 10001);

        Omega = 0.*t + Omega0;  % constant pulse
        Delta = 0.*t;           % zero detuning
        controls = [t',Delta',Omega'];

    else 
        disp(['The type of control can be: ''pi'', ''alphaRIO'', or ''alphadeltaRIO''']);
    end
end

