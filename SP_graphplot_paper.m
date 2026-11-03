% ==========================================================================
% SP_graphplot_paper: 論文の図 / the four figures of the ACNDC 2027 paper
%
%   result/figures/open_loop_loops.png                    (Fig. 1)
%   result/figures/relaxation_hysteresis.png              (Fig. 3)
%   result/figures/nominal_physical_states.png            (Fig. 4)
%   result/figures/nominal_commands_and_hj_selection.png  (Fig. 5)
% Compare with Documents/NS_DD_ACNDC_2027_figures/ (crops of the native run).
% P5_solution_check.m で使用
% ==========================================================================

	figureFolder = fullfile(cfg.resultFolder, 'figures');
	figureFiles = f_paper_figures(caseLQR, caseHJ, relaxCases, openCases, des, cfg, figureFolder);
	fprintf('\nPaper figures written to %s\n', figureFolder);
	for i_fg = 1:numel(figureFiles), fprintf('   %s\n', figureFiles{i_fg}); end
	clear i_fg
