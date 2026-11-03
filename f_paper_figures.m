function files = f_paper_figures(caseLQR, caseHJ, relaxCases, openCases, des, cfg, folder)
% ==================================================================
% files = f_paper_figures(caseLQR, caseHJ, relaxCases, openCases, des, cfg, folder)
%
% The four figures of the paper (Documents/NS_DD_ACNDC_2027.pdf), same file
% names as in its figures/ folder, written to <folder> as PNG (cfg.figureDPI)
% and, in MATLAB, as FIG:
%   open_loop_loops.png                    Fig. 1  theta vs sigma(u), h vs p
%   relaxation_hysteresis.png              Fig. 3  |h - h*| for each lambda_h
%   nominal_physical_states.png            Fig. 4  theta, v, p, h (LQR vs guarded HJ)
%   nominal_commands_and_hj_selection.png  Fig. 5  sigma(u), nu, ||y||, HJ selection
% (Fig. 2 of the paper is the TikZ block diagram of the sampled loop.)
% Contents, labels, colours and line styles follow the native figures from
% which the paper crops were taken; the figure sizes (inches) match the crops.
% Axes are placed explicitly (margins in inches) so that no label is clipped.
% ==================================================================
if ~exist(folder,'dir'), mkdir(folder); end
files = {};
colors = [0.35 0.35 0.35; 0.05 0.35 0.80];

% ---------------- Fig. 1: open-loop loops ----------------
k_sin = 0;
for k = 1:numel(openCases)
    if ~isempty(openCases{k}) && strcmp(openCases{k}.controller, 'sinusoidal'), k_sin = k; end
end
if k_sin > 0
    a = openCases{k_sin};
    sz = [6.6 1.75];
    fig = new_figure('Open-loop loops (Fig. 1)', cfg, sz);
    ax = grid_axes(fig, sz, 1, 2, [0.55 0.12 0.48 0.08 0.75 0]);
    plot(ax(1), a.u, a.x(:,1)); grid(ax(1), 'on'); xlim(ax(1), [0 1]);
    xlabel(ax(1), 'Applied command sigma(u)'); ylabel(ax(1), 'Angle theta');
    plot(ax(2), a.x(:,3), a.x(:,4)); grid(ax(2), 'on'); xlim(ax(2), [0 1]);
    xlabel(ax(2), 'Pressure'); ylabel(ax(2), 'Hysteresis h');
    files{end+1} = save_figure(fig, folder, 'open_loop_loops', cfg);
end

% ---------------- Fig. 3: relaxation of the hysteresis ----------------
if ~isempty(relaxCases)
    sz = [6.8 2.6];
    fig = new_figure('Relaxation sensitivity (Fig. 3)', cfg, sz);
    ax = grid_axes(fig, sz, 1, 1, [0.6 0.15 0.48 0.1 0 0]);
    hold(ax, 'on')
    for k = 1:numel(relaxCases)
        a = relaxCases{k};
        hS = 0; if isfield(a,'hStar'), hS = a.hStar; end
        plot(ax, a.t, abs(a.x(:,4) - hS), 'DisplayName', sprintf('lambda_h=%g', a.lambda_h));
    end
    grid(ax, 'on'); box(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, '|h-hStar|'); legend(ax, 'show');
    files{end+1} = save_figure(fig, folder, 'relaxation_hysteresis', cfg);
end

% ---------------- Fig. 4: nominal physical states ----------------
sz = [6.7 5.45];
fig = new_figure('Nominal physical states (Fig. 4)', cfg, sz);
ax = grid_axes(fig, sz, 2, 2, [0.6 0.15 0.45 0.3 0.75 0.8]);
labs = {'Bending angle (rad)','Bending velocity (rad/s)','Pressure (normalized)','Hysteresis h'};
cs = {caseLQR, caseHJ};
for j = 1:4
    hold(ax(j), 'on')
    for k = 1:2
        plot(ax(j), cs{k}.t, cs{k}.x(:,j), 'Color', colors(k,:), 'LineWidth', 1.2);
    end
    plot(ax(j), [0 caseLQR.t(end)], des.xStar(j)*[1 1], 'k:');
    grid(ax(j), 'on'); box(ax(j), 'on'); xlabel(ax(j), 'Time (s)'); ylabel(ax(j), labs{j}); title(ax(j), labs{j});
end
legend(ax(1), 'LQR', 'Guarded HJ', 'Equilibrium', 'Location', 'best');
files{end+1} = save_figure(fig, folder, 'nominal_physical_states', cfg);

% ---------------- Fig. 5: commands and HJ selection ----------------
sz = [6.7 5.4];
fig = new_figure('Nominal commands and HJ selection (Fig. 5)', cfg, sz);
ax = grid_axes(fig, sz, 4, 1, [0.75 0.15 0.45 0.1 0 0.3]);
hold(ax(1), 'on')
for k = 1:2, plot(ax(1), cs{k}.t, cs{k}.u, 'Color', colors(k,:)); end
plot(ax(1), [0 caseLQR.t(end)], cfg.uMax*[1 1], 'r:');
ylabel(ax(1), 'Applied sigma(u)'); grid(ax(1), 'on'); box(ax(1), 'on');
legend(ax(1), 'LQR', 'Guarded HJ', 'Upper bound');
hold(ax(2), 'on')
for k = 1:2, stairs(ax(2), cs{k}.t(1:end-1), cs{k}.nu, 'Color', colors(k,:)); end
ylabel(ax(2), 'Raw rate nu=du/dt'); grid(ax(2), 'on'); box(ax(2), 'on');
hold(ax(3), 'on')
for k = 1:2
    yk = (cs{k}.z - des.zStar')./des.scale';
    plot(ax(3), cs{k}.t, max(sqrt(sum(yk.^2,2)), 1e-12), 'Color', colors(k,:));
end
ylabel(ax(3), 'Scaled augmented error'); grid(ax(3), 'on'); box(ax(3), 'on');
stairs(ax(4), caseHJ.t(1:end-1), double(caseHJ.hjUsed)); ylim(ax(4), [-0.1 1.1]);
ylabel(ax(4), 'HJ selected'); xlabel(ax(4), 'Time (s)'); grid(ax(4), 'on');
files{end+1} = save_figure(fig, folder, 'nominal_commands_and_hj_selection', cfg);
end

% ------------------------------------------------------------------------------
function fig = new_figure(name, cfg, sizeInches)
fig = figure('Name', name, 'Visible', cfg.figureVisible, 'Color', 'w', ...
    'Units', 'inches', 'Position', [1 1 sizeInches]);
set(fig, 'PaperUnits', 'inches', 'PaperPositionMode', 'manual', ...
    'PaperPosition', [0 0 sizeInches], 'PaperSize', sizeInches);
end

function ax = grid_axes(fig, sizeInches, nr, nc, marg)
% marg = [left right bottom top hgap vgap] in inches; axes numbered row-wise
W = sizeInches(1); H = sizeInches(2);
w = (W - marg(1) - marg(2) - (nc-1)*marg(5))/nc;
h = (H - marg(3) - marg(4) - (nr-1)*marg(6))/nr;
ax = gobjects(nr*nc, 1);
for r = 1:nr
    for c = 1:nc
        k = (r-1)*nc + c;
        x0 = marg(1) + (c-1)*(w + marg(5));
        y0 = H - marg(4) - r*h - (r-1)*marg(6);
        ax(k) = axes('Parent', fig, 'Units', 'normalized', 'Position', [x0/W y0/H w/W h/H]);
    end
end
end

function file = save_figure(fig, folder, name, cfg)
file = fullfile(folder, [name '.png']);
if cfg.saveFigures
    print(fig, file, '-dpng', sprintf('-r%d', cfg.figureDPI));
    if ~exist('OCTAVE_VERSION','builtin') && exist('savefig','file')
        savefig(fig, fullfile(folder, [name '.fig']));
    end
end
end
