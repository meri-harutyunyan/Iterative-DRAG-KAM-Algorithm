function [RRobustnessDelta, RobustnessDelta2] = RobustnessDelta(DeltaDev, t, HH, maxN)
% -------------------------------------------------------------------------
% RobustnessDelta: Computes robustness of DRAG control vs detuning deviation
%
% Inputs:
%   DeltaDev - array of detuning deviations (fractional or scaled)
%   t        - time vector
%   HH       - Hamiltonian tensor (3x3x(maxN+1)xnTime)
%   maxN     - maximum DRAG order
%
% Outputs:
%   RRobustnessDelta - robustness values for 3-level system
%   RobustnessDelta2 - robustness values for 2-level reference
% -------------------------------------------------------------------------

dt = t(2) - t(1);
nTime = length(t);
psiT = [0; -1i; 0];  % target state

% Maximum control amplitude (used to scale detuning deviations)
maxOmega = max(abs(HH(2,1,2,:)))*2;
% ------------------- 3-Level Robustness -------------------- %
RRobustnessDelta = [];
for N = 1:maxN
    N
    RobustnessDelta = [];

    for dev = DeltaDev*maxOmega
        phi = [1; 0; 0];  % initial state
        for t_i = 1:nTime
            H = sum(HH(:,:,1:N+1,t_i), 3);  % sum Hamiltonian over orders
            H(2,2) = H(2,2) + dev;  % apply detuning deviation
            H(3,3) = H(3,3) + 2*dev;  % apply detuning deviation

            phi = expm(-1i * H * dt) * phi;
        end

        RobustnessDelta = [RobustnessDelta; abs(psiT' * phi)^2];
    end

    RRobustnessDelta = [RRobustnessDelta, RobustnessDelta];
end

% ------------------- 2-Level Reference -------------------- %
RobustnessDelta2 = [];
for dev = DeltaDev*maxOmega
    phi = [1; 0];
    for t_i = 1:nTime
        H = sum(HH(1:2,1:2,1:2,t_i), 3);
        H(2,2) = H(2,2) + dev;
        phi = expm(-1i * H * dt) * phi;
    end

    RobustnessDelta2 = [RobustnessDelta2; abs([0 -1i] * phi)^2];
end

end
