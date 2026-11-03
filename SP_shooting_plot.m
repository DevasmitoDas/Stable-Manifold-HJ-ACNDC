% ==========================================================================
% SP_shooting_plot: 計算結果をfigureで表示 / anchors on every (y_i, y_j) plane
% (dbFlag = 1). Circles: (y(0)) of each accepted trajectory.
% P5_solution_check.m で使用
% ==========================================================================

	nch = nchoosek(1:dim, 2);
	figure(2001); clf
	for i_sp = 1:size(nch,1)
		subplot(2, ceil(size(nch,1)/2), i_sp)
		plot(txpsumT(:,1+nch(i_sp,1)), txpsumT(:,1+nch(i_sp,2)), 'b.'); hold on
		plot(ini_sumT(:,nch(i_sp,1)), ini_sumT(:,nch(i_sp,2)), 'o'); hold off
		grid on
		xlabel(['y' num2str(nch(i_sp,1))]); ylabel(['y' num2str(nch(i_sp,2))]);
	end
	subplot(2, ceil(size(nch,1)/2), 1); legend('anchors (shooting lines)', 'y(0)')
	clear nch i_sp
