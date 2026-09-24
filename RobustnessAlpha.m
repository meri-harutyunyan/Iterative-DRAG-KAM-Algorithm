function [RRobustnessAlpha, RobustnessAlpha2] = RobustnessAlpha(AlphaDev, t, HH, maxN)
% -------------------------------------------------------------------------
% RobustnessAlpha: Computes robustness of DRAG control vs amplitude deviation
%
% Inputs:
%   AlphaDev   - array of amplitude deviations (e.g. -0.6:0.06:0.6)
%   t          - time vector
%   HH         - Hamiltonian tensor (3x3x(maxN+1)xnTime)
%   maxN       - maximum DRAG order
%
% Outputs:
%   RobustnessAlpha   - robustness values for 3-level system
%   RobustnessAlpha2  - robustness values for 2-level system
%
% -------------------------------------------------------------------------

dt = t(2) - t(1);
nTime = length(t);
psiT = [0; -1i; 0]; % target state

% ------------------- 3-Level Robustness -------------------- %
RRobustnessAlpha = [];
for N = 1:maxN
    N
    RobustnessAlpha = [];

    for dev = AlphaDev
        phi = [1; 0; 0]; % initial state
        
        for t_i = 1:nTime
           
            H = sum(HH(:,:,1:N+1,t_i), 3);  % Hamiltonian for order N

            % multiply coupling terms by (1 + dev)
            idx = [1 2; 2 1; 2 3; 3 2];
            H(sub2ind(size(H), idx(:,1), idx(:,2))) = ...
                (1 + dev) * H(sub2ind(size(H), idx(:,1), idx(:,2)));

            % time evolution
            phi = expm(-1i * H * dt) * phi;
        end

        % final fidelity for this deviation
        RobustnessAlpha = [RobustnessAlpha; abs(psiT' * phi)^2];
    end

    RRobustnessAlpha = [RRobustnessAlpha, RobustnessAlpha];
end

% ------------------- 2-Level Reference -------------------- %
RobustnessAlpha2 = [];
for dev = AlphaDev
    phi = [1; 0];
    
    for t_i = 1:nTime
        H = sum(HH(1:2,1:2,1:2,t_i), 3);

        idx = [1 2; 2 1];
        H(sub2ind(size(H), idx(:,1), idx(:,2))) = ...
            (1 + dev) * H(sub2ind(size(H), idx(:,1), idx(:,2)));

        phi = expm(-1i * H * dt) * phi;
    end
    
    RobustnessAlpha2 = [RobustnessAlpha2; abs([0 -1i] * phi)^2];
end

end
