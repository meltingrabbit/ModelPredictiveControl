function u = MPC(x, u0, A, B, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max)
	% モデル予測制御
	% u 		: 次のステップの出力					: NU
	% x 		: 状態変数								: NX
	% u0 		: 1ステップ前の入力						: NU
	% A 		: 離散化 状態方程式						: NX x NX
	% B 		: 離散化 状態方程式						: NX x NU
	% r 		: 目標状態（今は固定を想定している）	: NX
	% Hp		: 予測ホライズン						: 1
	% Hu		: 制御ホライズン						: 1
	% Q			: 状態量に対する重み行列 				: NX x NX
	% Qf		: 状態量に対する最終重み行列			: NX x NX
	% R			: 入力変化に対する重み行列				: NU x NU
	% x_min		: 制約条件								: NX
	% x_max		: 制約条件								: NX
	% u_min		: 制約条件								: NU
	% u_max		: 制約条件								: NU


	% ######################################################
	% MPC
	% ######################################################
	% cは飾り文字 mathcal
	% x(k+1) = A x(k) + B u(k)
	% y(k+1) = Cy x(k)				: 測定ベクトル
	% z(k+1) = Cz x(k)				: 出力ベクトル
	% r(k) 							: 目標値，参照軌道

	% p.90
	% cZ(k) 	: 1~Hpまでのz（ここでは =y =x）の予測ベクトル
	% cT(k) 	: 1~Hpまでのrの目標値ベクトル
	% ⊿U(k)	: 0~Huまでの⊿uの予測ベクトル

	% p.91
	% cZ(k) = Pis x(k) + Ups u(k-1) + Theta ⊿U(k)
	% cE(k) = cT(k) - Pis x(k) - Ups u(k-1)
	% ######################################################


	%% 初期化など
	persistent flag NX NU cQ cR;
	% flag 		: 初回呼び出し用						: 1
	% NX 		: xのサイズ								: 1
	% NU 		: uのサイズ								: 1
	% cは飾り文字 mathcal
	% cQ 		: 重み行列Q p.91 (3.3) 					: NX*Hp x NX*Hp
	% cR 		: 重み行列R p.91 (3.4) 					: NU*Hu x NU*Hu
	% Psi 		: Psi p.91 (3.5), p68					: NX*Hp x NX
	% Ups 		: Upsilon p.91 (3.5), p68				: NX*Hp x NU
	% Theta 	: Theta p.91 (3.5), p68					: NX*Hp x NU*Hu

	if isempty(flag)
		[flag, NX, NU, cQ, cR] ...
					=  MPC_init(x, u0, A, B, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max);
	end



	u = zeros(NU, 1);


end


function [flag, NX, NU, cQ, cR] = MPC_init(x, u0, A, B, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max)
	% 初回１度だけ呼ばれる準備関数


	[NX, NU] = size(B);
	cQ = BlkdiagNTime(Q, Hp-1);
	cQ = blkdiag(cQ, Qf);
	cR = BlkdiagNTime(R, Hu);

	Psi   = zeros(NX*Hp, NX);
	Ups   = zeros(NX*Hp, NU);
	Theta = zeros(NX*Hp, NU*Hu);

	temp_Apow    = eye(NX);
	temp_UpsPre  = zeros(NX, NU);
	for i = 1:Hp
		Ups(NX*(i-1)+1:NX*i, :) = temp_Apow * B + temp_UpsPre;
		temp_UpsPre = Ups(NX*(i-1)+1:NX*i, :);

		temp_Apow = temp_Apow * A;
		Psi(NX*(i-1)+1:NX*i, :) = temp_Apow;
	end



	flag = 1;
end

