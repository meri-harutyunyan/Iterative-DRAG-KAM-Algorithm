function Fid = DRAG(t,Delta1,Omega1,maxN,d0,kappa)
global HH
lambda=sqrt(2); % H_{23}/H_{12}
psiT=[0;-1i;0]; % target state for complete population transfer

dt=t(2)-t(1); 
nTime=length(t); % number of timesteps

% ------------------- Transformation W_1 -------------------- %
WW = zeros(3,3,nTime);
WW(1,2,:) = d0*kappa*Omega1/2;  WW(2,1,:) =  -d0*kappa*Omega1/2;
WW(2,3,:) = d0*lambda*Omega1/2; WW(3,2,:) =  -d0*lambda*Omega1/2;


% ------------------- Hamiltonian --------------------------- %
H0=[0 0 0; 0 0 0; 0 0 d0]; % H_0 without controls
HH = zeros(3,3,maxN+1,nTime); % HH(i,j,k,l) = H_{ij}^{(k-1)}(t_l)
HH(:,:,1,:) = repmat(H0,1,1,nTime);
HH(1,2,2,:) =  Omega1/2; HH(2,1,2,:) =  Omega1/2; HH(2,2,2,:) =  Delta1;
HH(2,3,2,:) =  lambda*Omega1/2; HH(3,2,2,:) =  lambda*Omega1/2;
HH(3,3,2,:) =  2*Delta1;
% ---------- Transformations and Dynamical Corrections ------ %
calTT = zeros(3, 3, maxN+1, nTime); % \cal{T}_N = T_1*T_2*...*T_N
calTT(:,:,1,:) = repmat(eye(3), 1, 1, nTime);
DDC = calTT; 

% ------------------- Fidelity at N=1 ----------------------- %
phi = [1;0;0]; % initial state
for t_i=1:nTime
    H = sum(HH(:,:,:,t_i),3); % H_1 = H_0 + V^(1)
    phi=expm(-1i*H*dt)*phi;
end
Fid = (abs(psiT'*phi))^2;

% ------------------- Main Algorithm ------------------------ %
for N=2:maxN
    N
    phi = [1;0;0]; % inital state
    [~, ~, dWW] = gradient(WW, 1, 1, dt); % derivative

    for t_i=1:nTime
        % ------------------- Transformations -------------------- %
        W = WW(:,:,t_i); % W_{N-1}
        dW = dWW(:,:,t_i); % derivative
        T = seriesT(N-1,W,maxN+1); % T_{N-1}
        T_d = conj(permute(T, [2 1 3])); % dagger
        dT= seriesdT(N-1,W,dW,maxN+1);
        calT = productSeries(calTT(:,:,:,t_i), T); % \cal{T}_{N-1} = T_1*T_2*...*T_{N-1}
        calT_d = conj(permute(calT, [2 1 3]));
        calTT(:,:,:,t_i) = calT;  % save \cal{T}_{N-1}

        % ------------------- DC_{N-1} calculation --------------- %
        DC = DDC(:,:,:,t_i);
        DC = productSeries(DC,T);
        DC = productSeries(T_d,DC);
        DC = DC + productSeries(T_d,dT);
        DDC(:,:,:,t_i) = DC; % save DC_{N-1}

        % ------------------- V^(N) calculation -------------------- %
        G = productSeries(calT_d,HH(:,:,:,t_i));
        G = productSeries(G,calT); % G_{N-1}
        OmegaN =  - G(2,1,N+1) + 1i*DC(2,1,N); 
        DeltaN =  - G(2,2,N+1) + 1i*DC(2,2,N) + G(1,1,N+1) - 1i*DC(1,1,N);
        HH(:,:,N+1,t_i) = [0       conj(OmegaN)          0;
                           OmegaN  DeltaN                lambda*conj(OmegaN);
                           0       lambda*OmegaN         2*DeltaN]; % V^(N)

        % ------------------- W_N calculation -------------------- %
        W_13 = G(1,3,N+1) - 1i*DC(1,3,N) + HH(1,3,N+1,t_i);
        W_23 = G(2,3,N+1) - 1i*DC(2,3,N) + HH(2,3,N+1,t_i);
        WW(:,:,t_i) = d0*[0 0 W_13; 0 0 W_23; -conj(W_13) -conj(W_23) 0]; % W_N

        % ------------------- Hamiltonian dynamics --------------- %
        H = sum(HH(:,:,:,t_i),3); % H_N
        phi=expm(-1i*H*dt)*phi; % Schrodinger's equation
        
    end
    Fid = [Fid (abs(psiT'*phi))^2]; % Fidelity at N
    
end

