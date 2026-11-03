% ==================================================================
% 非線形最適制御コントローラ算出プログラム / Nonlinear optimal controller design
%
% Saturation-aware stable-manifold optimal control of a hysteretic soft
% pneumatic actuator  (N. Sakamoto, D. Das, ACNDC 2027)
%
% - 安定多様体論を適用し, 漸化式を直接計算で解いて (y(t), mu(t)) を求める.
%   Stable-manifold trajectories by the integral iteration (P1-P4).
% - 価値関数 Vhat(y) を多様体と HJB 残差から求める (P5).
% - ガード付き HJ / LQR 切替の性能を確認し, 論文の図を作成する (P5).
%
% Run from the repository root:   >> DirectApp
% Outputs (result/):  <date>_..._data.mat (after P4), <date>_..._results.mat,
%   figures/open_loop_loops.png, relaxation_hysteresis.png,
%   nominal_physical_states.png, nominal_commands_and_hj_selection.png,
%   numerical_summary.txt, *.csv
%
% ※）評価関数 J には1/2を付けない形でQ,Rを与える. (Q, R without the factor 1/2)
% ==================================================================

clear all %#ok<CLALL>
close all
clc

% === パスの追加 / paths ===
rootFolder = fileparts(mfilename('fullpath'));
if isempty(rootFolder), rootFolder = pwd; end
cd(rootFolder)
addpath(rootFolder)
addpath(fullfile(rootFolder,'inputData'))
addpath(fullfile(rootFolder,'polybase'))
addpath(fullfile(rootFolder,'VariousFunction'))
addpath(fullfile(rootFolder,'SP_iterData'))

% == problem to be solved == 解きたい問題のフォルダを指定 == %
addpath(fullfile(rootFolder,'inputData','5D_Pneumatic_nu_sat_symb'))

% === run options ===
askFlag = 0;                   % 1: interactive prompts as in the reference DirectApp
tmp_FileName = 'ACNDC2027';    % tag of the saved MAT files
if askFlag == 1
	tmp_FileName = input('save MAT-file as ...  ','s');
	while isempty(tmp_FileName)
		tmp_FileName = input('save MAT-file as ...  ','s');
	end
end

% === シミュレーション結果を格納するフォルダ作成 / result folder ===
if ~exist(fullfile(rootFolder,'result'),'dir')
	mkdir(fullfile(rootFolder,'result'));
end

P1_data_setting      % プログラム1 : 共通のパラメータの設定 / local design
P2_hamiltonian_check % プログラム2 : 収束の確認 / convergence along eigen-directions
P3_decide_xi         % プログラム3 : 初期値ξ / seeds xi

if askFlag == 1
	fprintf('\nPress any key to continue\n')
	pause
end

tic; % 繰り返し計算の計算時間計測開始
P4_iterative_calcu   % プログラム4 : ξに対して(y(t), mu(t))を求める繰り返し演算
loopend1 = toc;

% === 計算時間の表示 ===
nhour1 = floor(loopend1/60/60);
nmin1 = floor((loopend1 - nhour1*3600)/60);
nsec1 = loopend1 - nhour1*3600 - nmin1*60;
fprintf('%s%d%s%02d%s%04.1f%s\n', 'computation time for the stable manifold : ', ...
	nhour1,'h ',nmin1,'min ',nsec1,'sec');

% === 計算結果の保存1 / save data up to now (reference file naming) ===
Today = clock;
FileName1 = sprintf('%d_%d_%d_%d%d_dim=%d_k=%d_radi=%g_%s_data.mat', ...
	Today(1),Today(2),Today(3),Today(4),Today(5),dim,knum,radi,tmp_FileName);
save(fullfile(rootFolder,'result',FileName1),'-v7');

P5_solution_check    % プログラム5 : 価値関数, 閉ループ, 論文の図 / value fit, loops, figures

% === 計算結果の保存2 / save everything ===
FileName2 = strrep(FileName1,'_data.mat','_results.mat');
save(fullfile(rootFolder,'result',FileName2),'-v7');
fprintf('\nResults saved in %s\n', fullfile(rootFolder,'result'));
