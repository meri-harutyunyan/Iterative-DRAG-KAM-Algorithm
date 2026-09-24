clear all
format long

global HH

%-----Controlable variables-----%
kappa=0.5; % beta/alpha
n=2; % degree of hyperGaussian (Omega)
nTime = 1000; % number of timesteps after conversion to a hyperGaussian
nsigma = 5; % Duration/(2*sigma)
maxN = 8; % Number of iterations
d2=-1; % H_33, 1 or -1 due to normalization
ssig = [0.5 0.75 1.5 2]*2*pi/abs(d2); % Sigma of the hyperGaussian Omega

TypeOfControl = 'pi'; % Can be 'pi', 'alphaRIO', 'alphadeltaRIO'
PlotRobustnessAlpha = 0;
PlotRobustnessDelta = 0;

% Getting the controls with constant Omega
controls = getControls(TypeOfControl);
t = controls(:,1);
Delta = controls(:,2);
Omega = controls(:,3); 
% ! Or simply read any controls as t, Delta, Omega !


FFid = [];
for sig=ssig
    [tau Delta1 Omega1] = convertControls(t,Delta,Omega,sig,n,nTime,nsigma); % Converting the controls so that Omega is a hyperGaussian
    Fid = DRAG(tau,Delta1,Omega1,maxN,d2,kappa);
    FFid = [FFid;Fid];
end

%%

N_list=[1:1:maxN];

figure('PaperUnits','centimeters','PaperPosition',[1,1,20,14]);
plot(N_list, log10(1-FFid'), '-'); 
ylim([-10, -1])
xlim([1 maxN])
ylabel("Infidelity",'Interpreter', 'latex')
xlabel('$N$, number of iterations','Interpreter', 'latex')
yticks(linspace(-14, 0, 8));
xticks(N_list);
set(gca, 'FontSize', 20);
legendInfo = arrayfun(@(x) ['$\sigma = ', sprintf('%.1f', x), '\sigma_0$'], ssig/2/pi, 'UniformOutput', false);
legend([legendInfo], 'Interpreter', 'latex', 'Location', 'southwest');
grid on

%%

if PlotRobustnessAlpha
    % ------------------- Parameters ------------------- %
    AlphaDev = -0.6:0.01:0.6;
    
    % ------------------- Compute Robustness ------------------- %
    [RRobustnessAlpha, RobustnessAlpha2] = RobustnessAlpha(AlphaDev, tau, HH, maxN);
    
    % ------------------- Plot ------------------- %
    figure('PaperUnits','centimeters','PaperPosition',[1,1,10,7]);
    axes('Position',[0.14,0.20,0.72,0.7])
    
    plot(AlphaDev, log10(1 - RRobustnessAlpha))
    hold on
    plot(AlphaDev, log10(1 - RobustnessAlpha2), 'k--')
    
    legendInfo = cellstr(num2str((1:maxN)', '$N= %d$'));
    legend([legendInfo; "2 levels"], 'Location', 'southeast', 'Interpreter', 'latex');
    
    xlabel('$\alpha$, Amplitude deviation', 'Interpreter', 'latex')
    ylabel('Infidelity', 'Interpreter', 'latex')
    xlim([-0.42, 0.42])
    ylim([-4.6, -0.5])
    yticks([-4 -3 -2 -1])
    set(gca, 'FontSize', 10)
    grid on
end





%% --- Plot Robustness vs Detuning Deviation (Delta) --- %%
if PlotRobustnessDelta
    % ------------------- Parameters ------------------- %
    DeltaDev = -0.6:0.01:0.6;
    
    % ------------------- Compute Robustness ------------------- %
    [RRobustnessDelta, RobustnessDelta2] = RobustnessDelta(DeltaDev, tau, HH, maxN);
    
    % ------------------- Plot ------------------- %
    figure('PaperUnits','centimeters','PaperPosition',[1,1,10,7]);
    axes('Position',[0.14,0.20,0.72,0.7])
    
    plot(DeltaDev, log10(1 - RRobustnessDelta))
    hold on
    plot(DeltaDev, log10(1 - RobustnessDelta2), 'k--')
    
    legendInfo = cellstr(num2str((1:maxN)', '$N= %d$'));
    legend([legendInfo; "2 levels"], 'Location', 'southeast', 'Interpreter', 'latex');
    
    xlabel('$\delta/\Omega_{\max}$, Detuning deviation', 'Interpreter', 'latex')
    ylabel('Infidelity', 'Interpreter', 'latex')
    xlim([-0.42, 0.42])
    ylim([-4.6, -0.5])
    yticks([-4 -3 -2 -1])
    set(gca, 'FontSize', 20)
    grid on
end
