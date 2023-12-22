%% Sanity check 2: Does the Eq. (vi) recover the velocity source condition?
% (nondimensionalized form is Eq. (vii))
close all; clear all; clc ; 
set(groot,'DefaultAxesFontSize', 15); set(groot,'DefaultLineLineWidth', 2); set(groot,'DefaultAxesLineWidth', 2)

gamma = 1/2; %gamma = a/b
R = linspace(0,1/gamma); % R = r/a
maxn = 40; %maximum term in summation 
summation = 0; %initialize summation
for n = 1:maxn
    summation = summation + 1/(BessDerivZerosBisect2(0,n))*besselj(1,BessDerivZerosBisect2(0,n)*gamma)...
        /(besselj(0,BessDerivZerosBisect2(0,n)))^2*besselj(0,besselj(0,BessDerivZerosBisect2(0,n))*R*gamma);
end
LHS = 2*gamma*summation; %LHS of Eq. (vii)
RHS = 1.*(R<=1.0); %RHS of Eq. (vii)
plot(R,LHS); hold on; plot(R,RHS);
xlim([min(R),max(R)]); ylim([0,1]); xlabel('$R \equiv r/a$','Interpreter','latex');
ptitle = sprintf('First %0.5g terms of Eq. (vii)', maxn); title(ptitle);
legend(['$2 \gamma \sum_{n=1}^\infty \frac{1}{\alpha_{0n}^\prime}' ...
    ' \frac{J_1(\alpha_{0n}^\prime\gamma)}{J_0^2(\alpha^\prime_{0n}) } J_0(\alpha^\prime_{0n} \gamma R)$'],...
    '$1, R \in [0,1];\quad 0, R \in (1, 1/\gamma]$','Interpreter','latex')

% Add these lines to save the figure as a transparent SVG file using export_fig
filename = 'first_40_terms.svg'; % Set your desired filename
export_fig(filename, '-transparent');
