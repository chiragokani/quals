%% Sanity check 2: Does the Eq. (vi) recover the velocity source condition?
% (nondimensionalized form is Eq. (vii))
close all; clear all; clc ; 
set(groot,'DefaultAxesFontSize', 15); set(groot,'DefaultLineLineWidth', 2); set(groot,'DefaultAxesLineWidth', 2)

gamma = 1/2; %gamma = a/b
R = linspace(0,1/gamma,500); % R = r/a
maxn = 40; %maximum term in summation 
alphap = alpha0np(maxn);
alphap(1) = 1e-9; % Make very small instead of 0 for numerical stability

summation = 0; %initialize summation
for n = 1:maxn
    summation = summation + 1/alphap(n)*besselj(1,alphap(n)*gamma)...
        /(besselj(0,alphap(n)))^2*besselj(0,alphap(n)*R*gamma);
end
LHS = 2*gamma*summation; %LHS of Eq. (vii)
RHS = 1.*(R<=1.0); %RHS of Eq. (vii)
plot(R,LHS); hold on; plot(R,RHS);
xlim([min(R),max(R)]); ylim([-0.1,1.1]); xlabel('$R \equiv r/a$','Interpreter','latex');
ptitle = sprintf('First %0.5g terms of Eq. (vii)', maxn); title(ptitle);
legend(['$2 \gamma \sum_{n=1}^\infty \frac{1}{\alpha_{0n}^\prime}' ...
    ' \frac{J_1(\alpha_{0n}^\prime\gamma)}{J_0^2(\alpha^\prime_{0n}) } J_0(\alpha^\prime_{0n} \gamma R)$'],...
    '$1, R \in [0,1];\quad 0, R \in (1, 1/\gamma]$','Interpreter','latex')
legend('location','east')
% Add these lines to save the figure as a transparent SVG file using export_fig
filename = '    first_40_terms.png'; % Set your desired filename
export_fig(filename, '-transparent');

%% Calculate alphap = [alpha01', alpha02', ..., alpha0N'] (by Jack Hallveld)
function alphap = alpha0np(N)
% alpha0n' is the nth root of J0'(x) = 0.
alphap = zeros(1,N);
step = 2;
xmin = 0.1;
J0p = @(x) -besselj(1,x); % J0'(x) = -J1(x)
n = 2; % Skip n=1 because already know alphap(1) = 0
while n <= N 
    % Look for the first zero crossing over the vector x
    x = linspace(xmin, xmin+step, 1000);
    y = J0p(x);
    ii = 1;
    while ii < length(x)
        if y(ii) * y(ii+1) <= 0
            % Found zero crossing; use fzero to zone in on the root
            alphap(n) = fzero(J0p, [x(ii), x(ii+1)]);
            xmin = x(ii+1);
            n = n+1;
            break
        end
        ii = ii+1;
    end
    if ii == length(x)
        % Went through the entire vector x without finding anything, so
        % need to look further out
        xmin = xmin+step;
    end
end
end
