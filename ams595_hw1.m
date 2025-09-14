
clear; clc; close all;
rng(0,'twister');   % Reproducibility

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose a maximum total number of points
Nmax = 1e6;                        
% sample checkpoints
idxToShow = unique(round(logspace(2, log10(Nmax), 300))); 

% Generate all points up-front for speed, then analyze progressively
x = rand(Nmax,1);
y = rand(Nmax,1);
inside = (x.^2 + y.^2) <= 1;

% Running estimates of pi as N grows from 1 to Nmax
cumInside = cumsum(inside);
Nvec       = (1:Nmax)';
piRunning  = 4 * (cumInside ./ Nvec);

% Values at sampled checkpoints for plotting
pi_at_idx  = piRunning(idxToShow);
err_at_idx = abs(pi_at_idx - pi);

% Plot: running estimate and deviation from true pi
figure('Name','Task 1: Running estimate and deviation','Color','w');
tiledlayout(2,1,"TileSpacing","compact","Padding","compact");

% Top: estimated pi vs N (with true pi reference)
nexttile;
semilogx(idxToShow, pi_at_idx, 'b-', 'LineWidth', 1.0); hold on;
yline(pi, 'r--', 'LineWidth', 1.2);
xlabel('Number of points N');
ylabel('Estimated \pi');
title('Running Monte Carlo estimate of \pi (for-loop, fixed N)');
legend('Estimate','True \pi','Location','best');

% Bottom: absolute deviation vs N
nexttile;
loglog(idxToShow, max(err_at_idx, eps), 'k-', 'LineWidth', 1.0);
xlabel('Number of points N (log scale)');
ylabel('|Estimate - \pi| (log scale)');
title('Absolute deviation vs N');

% Timing study: precision vs computational cost
% Evaluate a set of different total point counts and time each run.
Ngrid = round(logspace(3,7,14));   % from 1e3 to 1e7
time_per_N = zeros(size(Ngrid));
err_per_N  = zeros(size(Ngrid));
for i = 1:numel(Ngrid)
    Ni = Ngrid(i);
    tic;
    xi = rand(Ni,1);
    yi = rand(Ni,1);
    countInside = sum((xi.^2 + yi.^2) <= 1);
    pi_est = 4 * countInside / Ni;
    time_per_N(i) = toc;
    err_per_N(i)  = abs(pi_est - pi);   % using true pi only for evaluation/plotting
end

saveas(gcf, 'task1_running_estimate.png');

figure('Name','Task 1: Precision vs cost','Color','w');
loglog(time_per_N, max(err_per_N, eps), 'o-','LineWidth',1.0,'MarkerSize',6);
grid on;
xlabel('Execution time (s, log scale)');
ylabel('Absolute error |estimate - \pi| (log scale)');
title('Precision vs computational cost (for-loop, fixed N)');

saveas(gcf, 'task1_precision_vs_cost.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The stopping rule uses an uncertainty-based criterion:
% Require relative CI half-width (95% by default) <= 0.5*10^(-k),
% where k is the requested number of significant figures.

digitsList = 1:4;      % Selected precisions to achieve
alpha = 0.05;          % Two-sided significance level for CI (~95% coverage)
z = sqrt(2) * erfcinv(alpha);  % Normal critical value

batchSize = 5000;      % Generate points in batches for efficiency
fprintf('\nTask 2: While-loop precision run (no true pi used in stopping rule)\n');
fprintf('alpha = %.3f (two-sided), z = %.5f\n', alpha, z);
results = struct('digits',[],'pi_est',[],'N',[],'batches',[]);

for k = digitsList
    targetRel = 0.5 * 10^(-k);   % Relative half-width threshold for k significant figures

    N = 0;                      % Total points so far
    nInside = 0;                % Total points inside quarter circle
    batches = 0;                % While-loop iterations (batches)
    reached = false;

    % While-loop continues until the relative CI half-width criterion is met
    while ~reached
        % Generate a batch of random points
        xb = rand(batchSize,1);
        yb = rand(batchSize,1);
        nInside = nInside + sum((xb.^2 + yb.^2) <= 1);
        N = N + batchSize;
        batches = batches + 1;

        % Current estimate of pi and its standard error via binomial model
        p_hat = nInside / N;                                    % quarter-circle area
        pi_hat = 4 * p_hat;                                     % estimate of pi
        var_p = max(p_hat * (1 - p_hat), eps) / N;              % guard against zero
        se_pi = 4 * sqrt(var_p);                                % std error of pi_hat

        % Two-sided CI half-width and relative half-width (no true pi here)
        halfWidth = z * se_pi;
        relHalfWidth = halfWidth / max(pi_hat, eps);

        % Check precision criterion for k significant figures
        reached = (relHalfWidth <= targetRel);
    end

    % Record results
    results(end+1) = struct('digits',k,'pi_est',pi_hat,'N',N,'batches',batches); %#ok<SAGROW>
    fprintf('  k = %d sig figs: pi_hat = %.10f, N = %d points, batches = %d, relHalfWidth = %.3e\n', ...
            k, pi_hat, N, batches, relHalfWidth);
end
results = results(2:end); % remove the empty first struct element

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example calls:
[pi_k2, N_k2, it_k2] = mcpi_precision_plot(2);                         % simple call with default settings
[pi_k3, N_k3, it_k3] = mcpi_precision_plot(3,'BatchSize',10000);       % larger batches
[pi_k4, N_k4, it_k4] = mcpi_precision_plot(4,'Alpha',0.01,'RNGSeed',1);% 99% CI target

% The function required for Task 3 is defined below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [pi_hat, N, batches] = mcpi_precision_plot(nSigFigs, varargin)
% mcpi_precision_plot
% Estimate pi via Monte Carlo until a target significant-figure precision is met,
% using a while-loop and an uncertainty-based stopping rule that does not use true pi.
% A live plot shows random points as they are generated and saves the final plot.
%
% Usage:
%   [pi_hat, N, batches] = mcpi_precision_plot(nSigFigs)
%   [pi_hat, N, batches] = mcpi_precision_plot(nSigFigs, 'Name', Value, ...)
%
% Inputs:
%   nSigFigs  - integer, requested significant figures (e.g., 2, 3, 4)
%
% Name-Value options:
%   'Alpha'          - two-sided significance level for CI (default 0.05 for ~95% CI)
%   'BatchSize'      - number of new points per while-loop iteration (default 5000)
%   'MaxPoints'      - safety cap on total points to avoid endless runs (default 2e7)
%   'MaxPlotPoints'  - cap on total points stored for plotting (default 2e5)
%   'RNGSeed'        - set a seed for reproducibility (default [])
%   'MarkerSize'     - size of scatter markers (default 8)
%   'ShowPlot'       - true/false to show live plot (default true)
%   'SavePlot'       - true/false to save final plot (default true)
%   'PlotFormat'     - format for saving plot: 'png', 'fig', 'pdf', 'eps' (default 'png')
%   'PlotResolution' - DPI resolution for raster formats (default 300)
%   'OutputDir'      - directory to save plots (default current directory)
%
% Outputs:
%   pi_hat  - final estimate of pi
%   N       - total number of random points used
%   batches - number of while-loop iterations (batches)

    %---------------------- Parse inputs ----------------------%
    p = inputParser;
    p.addRequired('nSigFigs', @(v) validateattributes(v,{'numeric'},{'scalar','integer','>=',1}));
    p.addParameter('Alpha', 0.05, @(v) isnumeric(v) && isscalar(v) && v>0 && v<1);
    p.addParameter('BatchSize', 5000, @(v) isnumeric(v) && isscalar(v) && v>=1);
    p.addParameter('MaxPoints', 2e7, @(v) isnumeric(v) && isscalar(v) && v>=1);
    p.addParameter('MaxPlotPoints', 2e5, @(v) isnumeric(v) && isscalar(v) && v>=1);
    p.addParameter('RNGSeed', [], @(v) isempty(v) || (isscalar(v) && isnumeric(v)));
    p.addParameter('MarkerSize', 8, @(v) isnumeric(v) && isscalar(v) && v>0);
    p.addParameter('ShowPlot', true, @(v) islogical(v) && isscalar(v));
    p.addParameter('SavePlot', true, @(v) islogical(v) && isscalar(v));
    p.addParameter('PlotFormat', 'png', @(v) ismember(v, {'png', 'fig', 'pdf', 'eps', 'jpg', 'tiff'}));
    p.addParameter('PlotResolution', 300, @(v) isnumeric(v) && isscalar(v) && v>0);
    p.addParameter('OutputDir', '.', @(v) ischar(v) || isstring(v));
    p.parse(nSigFigs, varargin{:});
    S = p.Results;

    if ~isempty(S.RNGSeed)
        rng(S.RNGSeed, 'twister');
    end

    % Create output directory if it doesn't exist
    if S.SavePlot && ~isfolder(S.OutputDir)
        mkdir(S.OutputDir);
    end

    %----------------- Precision/CI configuration --------------%
    targetRel = 0.5 * 10^(-nSigFigs);     % relative CI half-width threshold
    z = sqrt(2) * erfcinv(S.Alpha);       % normal quantile for two-sided alpha

    %---------------------- Figure setup -----------------------%
    if S.ShowPlot
        fig = figure('Name',sprintf('Monte Carlo \x03c0 to %d significant figures', nSigFigs), ...
               'Color','w', 'Position', [100, 100, 800, 700]);
        hold on; axis equal; box on;
        xlim([0 1]); ylim([0 1]);
        xlabel('x'); ylabel('y');
        title(sprintf('Monte Carlo estimation of \\pi (target: %d significant figures)', nSigFigs));

        % Draw quarter-circle boundary of radius 1
        th = linspace(0, pi/2, 200);
        plot(cos(th), sin(th), 'k-', 'LineWidth', 1.2);

        % Create scatter objects for points inside and outside
        hIn  = scatter(nan, nan, S.MarkerSize, [0.1 0.6 0.1], 'filled', 'DisplayName','Inside');
        hOut = scatter(nan, nan, S.MarkerSize, [0.85 0.2 0.2], 'filled', 'DisplayName','Outside');
        legend('Quarter circle','Inside','Outside','Location','southoutside','NumColumns',3);

        % Pre-allocate arrays for plotted points
        xin = []; yin = [];
        xout = []; yout = [];
    end

    %---------------- While-loop simulation --------------------%
    N = 0;               % total points
    nInside = 0;         % count inside quarter circle
    batches = 0;         % while-loop iterations
    pi_hat = NaN;        % final estimate

    reached = false;
    while ~reached
        % Safety check: stop if exceeding MaxPoints
        if N >= S.MaxPoints
            warning('Reached MaxPoints = %g before meeting the precision target.', S.MaxPoints);
            break;
        end

        % Draw a batch of random points
        xb = rand(S.BatchSize,1);
        yb = rand(S.BatchSize,1);
        in = (xb.^2 + yb.^2) <= 1;

        % Update counters
        nInside = nInside + sum(in);
        N = N + S.BatchSize;
        batches = batches + 1;

        % Update running estimate and uncertainty
        p_hat = nInside / N;
        pi_hat = 4 * p_hat;

        var_p = max(p_hat * (1 - p_hat), eps) / N;  % protect against zero
        se_pi = 4 * sqrt(var_p);
        halfWidth = z * se_pi;
        relHalfWidth = halfWidth / max(pi_hat, eps);

        % Update plot
        if S.ShowPlot
            % Append recent batch to plotting buffers
            xin  = [xin;  xb(in)];
            yin  = [yin;  yb(in)];
            xout = [xout; xb(~in)];
            yout = [yout; yb(~in)];

            % Cap for memory/plotting speed
            cap = S.MaxPlotPoints;
            if numel(xin) > cap,  xin  = xin(end-cap+1:end);  yin  = yin(end-cap+1:end);  end
            if numel(xout)> cap,  xout = xout(end-cap+1:end); yout = yout(end-cap+1:end); end

            % Push to graphics
            set(hIn,  'XData', xin,  'YData', yin);
            set(hOut, 'XData', xout, 'YData', yout);

            % Display a status annotation (updates each batch)
            statusStr = sprintf('N = %d, \\pi \\approx %.6f, rel half-width \\approx %.2e', ...
                                 N, pi_hat, relHalfWidth);
            if ~exist('hText','var') || ~isgraphics(hText)
                hText = text(0.02, 0.97, statusStr, 'Units','normalized', ...
                             'VerticalAlignment','top','FontName','monospace', ...
                             'FontSize',10,'BackgroundColor',[1 1 1 0.6]);
            else
                set(hText,'String',statusStr);
            end
            drawnow limitrate;
        end

        % Check precision target
        reached = (relHalfWidth <= targetRel);
    end

    %----------------------- Final output ----------------------%
    % Print to command window with the requested significant-figure formatting
    fmt = sprintf('%%.%dg', nSigFigs);          
    pi_str = sprintf(fmt, pi_hat);
    msg = sprintf('Monte Carlo pi to %d significant figures: %s (N = %d, batches = %d)', ...
                  nSigFigs, pi_str, N, batches);
    disp(msg);

    % Print on the plot
    if S.ShowPlot
        txt = sprintf('\\pi \\approx %s (N = %d)', pi_str, N);
        text(0.50, 0.10, txt, 'Units','normalized', 'FontWeight','bold', ...
             'HorizontalAlignment','center','FontSize',12, ...
             'BackgroundColor',[1 1 1 0.7], 'EdgeColor',[0.3 0.3 0.3]);
        
        % Save the plot 
        if S.SavePlot
            filename = fullfile(S.OutputDir, sprintf('pi_%d_sigfigs_N%d', nSigFigs, N));
            
            switch lower(S.PlotFormat)
                case 'png'
                    print(fig, filename, '-dpng', sprintf('-r%d', S.PlotResolution));
                case 'jpg'
                    print(fig, filename, '-djpeg', sprintf('-r%d', S.PlotResolution));
                case 'tiff'
                    print(fig, filename, '-dtiff', sprintf('-r%d', S.PlotResolution));
                case 'pdf'
                    print(fig, filename, '-dpdf', '-bestfit');
                case 'eps'
                    print(fig, filename, '-depsc', sprintf('-r%d', S.PlotResolution));
                case 'fig'
                    saveas(fig, filename, 'fig');
            end
            
            fprintf('Plot saved as: %s.%s\n', filename, S.PlotFormat);
        end
    end
end


