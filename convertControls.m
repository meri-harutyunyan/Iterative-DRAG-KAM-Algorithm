function [tau, Delta_r, Omega_r] = convertControls(t, Delta, Omega, sigma, n, Nsteps, nsigma)
%   CONVERT_CONTROLS uses time reparametrization technique
%   in order to convert any pulse to a hyperGaussian shape.
%
%   This function implements the time reparametrization procedure described in
%   the DRAG-KAM framework.
%   A new time variable τ is introduced via t = g(τ), ensuring that the
%   reparametrized Rabi frequency Ω_r(τ) vanishes at the pulse boundaries
%   while the total pulse area is preserved.
%
%   Theory:
%       Ω_r(τ) = g'(τ) · Ω(g(τ)) = Ω_r_0*exp(-(τ/σ)^n)
%       Δ_r(τ) = g'(τ) · Δ(g(τ))
%
%       The mapping g(τ) is determined from the area condition:
%           ∫_g_i^g Ω(t) dt = Ω_r_0 ∫_τ_i^τ exp(-(s/σ)^n) ds
%
%   Inputs:
%       t        - time array
%       Delta    - detuning Δ(t)
%       Omega    - Rabi frequency Ω(t)
%       sigma    - σ Gaussian width
%       n        - hyper-Gaussian order
%       Nsteps   - number of τ discretization points 
%       nsigma   - τ_max - τ_min/2σ 
%
%   Outputs:
%       tau     - rescaled time grid
%       Delta_r - rescaled detuning Δ_r(τ)
%       Omega_r - rescaled Rabi frequency Ω_r(τ)
%
%   Author: Meri Harutyunyan

% ---------------- Time grid in τ ----------------
dtau = 2 * sigma * nsigma / Nsteps;
tau  = (-sigma*nsigma : dtau : sigma*nsigma - dtau) + dtau/2;

% ---------------- Interpolation splines ----------------
t = t - t(1);
spline_Omega = spline(t, Omega);
spline_Delta = spline(t, Delta);

% ---------------- Hyper-Gaussian amplitude ----------------
Omega_r0 = trapz(t, Omega) * n / (2 * sigma * gamma(1/n));

% ---------------- Solve for g(τ) ----------------
g = zeros(size(tau));
Area_r = cumtrapz(tau, Omega_r0 * exp(-(tau ./ sigma).^n)); % Partial rescaled area 
for i = 2:length(tau)
    area_diff   = @(x) trapz(0:x/1000:x, ppval(spline_Omega, 0:x/1000:x)) - Area_r(i);
    g(i)        = fzero(area_diff, 3);
end

% ---------------- Rescaled controls ----------------
gdot    = gradient(g, dtau);
Omega_r = gdot .* ppval(spline_Omega, g);
Delta_r = gdot .* ppval(spline_Delta, g);
end