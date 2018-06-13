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
	% cU(k)		: 0~Huまでのuの予測ベクトル
	% ⊿U(k)	: 0~Huまでの⊿uの予測ベクトル

	% p.91
	% cZ(k) = Psi x(k) + Ups u(k-1) + Theta ⊿U(k)
	% cE(k) = cT(k) - Psi x(k) - Ups u(k-1)

	% Ωθ ≦ ω: 制約一次不等式作成
	% ######################################################


	%% 初期化など
	persistent flag NX NU cQ cR Psi Ups Theta F f G g BlockF;
	% flag 		: 初回呼び出し用						: 1
	% NX 		: xのサイズ								: 1
	% NU 		: uのサイズ								: 1
	% cは飾り文字 mathcal
	% cQ 		: 重み行列Q p.91 (3.3) 					: NX*Hp x NX*Hp
	% cR 		: 重み行列R p.91 (3.4) 					: NU*Hu x NU*Hu
	% Psi 		: Psi p.91 (3.5), p.68					: NX*Hp x NX
	% Ups 		: Upsilon p.91 (3.5), p.68				: NX*Hp x NU
	% Theta 	: Theta p.91 (3.5), p.68				: NX*Hp x NU*Hu
	% F 	: cUの制約条件行列 p.99 (3.35) のfを除いたもの 		: 2*NU*Hu x NU*Hu 制約条件は上限下限なので2*NU*Hu個
	% f 	: cUの制約条件ベクトル p.99 (3.35) のf 				: 2*NU*Hu x 1
	% G 	: cZの制約条件行列 p.99 (3.36) のgを除いたもの 		: 2*NX*Hp x NX*Hp 制約条件は上限下限なので2*NX*Hp個
	% g 	: cZの制約条件ベクトル p.99 (3.36) のg 				: 2*NX*Hp x 1
	% BlockF: Fの累乗したやつ p.99の下から3行目 				: 2*NU*Hu x NU*Hu
	if isempty(flag)
		[flag, NX, NU, cQ, cR, Psi, Ups, Theta, F, f, G, g, BlockF] ...
					=  MPC_init(A, B, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max);
	end



	%% 制約一次不等式生成
	% Ωθ ≦ ω: 制約一次不等式作成
	% Omega 	: Ω p.100 (3.41)の左辺行列				:
	% omega 	: ω p.100 (3.41)の右辺ベクトル			:
	Omega = vertcat(BlockF, G*Theta);
	omega = vertcat(-BlockF(:,1:NU)*u0-f, -G*(Psi*x + Ups*u0)-g);


	%% 参照軌道からの残渣，追従誤差
	cT = repmat(r, Hp, 1);
	cE = cT - Psi*x - Ups*u0;					% ε．p.91 (3.6)

	cG = 2 * Theta.' * cQ * cE;					% G p.91 (3.11)
	cH = 2 * (Theta.' * cQ * Theta + cR);		% 2*H p.91 (3.12)

	% これつけないと，quadprogで
	% 警告: ヘッシアンが対称ではありません。H=(H+H')/2 をリセットしています。
	% って出てくる．
	cH = round(cH, 10, 'significant');

	%% 最適化．二次計画法
	opts = optimoptions(@quadprog,'Display','off');
	% opts = optimoptions(@quadprog,'Display','iter-detailed');
	% opts = optimoptions(@quadprog);
	dU = quadprog(cH, -cG.', Omega, omega, [], [], [], [], [], opts);
	du = dU(1:NU,1);

	u = u0 + du;
	% u = zeros(NU, 1);
end


%% 初期化関数
% 初回だけ呼ばれるため，変化する変数は引数に取れない
function [flag, NX, NU, cQ, cR, Psi, Ups, Theta, F, f, G, g, BlockF] ...
		 = MPC_init(A, B, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max)
	% 初回１度だけ呼ばれる準備関数


	[NX, NU] = size(B);


	%% 重み行列 p.91
	cQ = BlkdiagNTime(Q, Hp-1);
	cQ = blkdiag(cQ, Qf);
	cR = BlkdiagNTime(R, Hu);


	%% 予測行列 p.91
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
	for i = 1:Hu
		Theta( NX*(i-1)+1:NX*Hp , NU*(i-1)+1:NU*i ) = Ups( 1:NX*(Hp-i+1) , : );
	end


	%% 制約条件 p.99
	F = zeros(2*NU*Hu,NU*Hu);
	f = zeros(2*NU*Hu,1);
	G = zeros(2*NX*Hp,NX*Hp);
	g = zeros(2*NX*Hp,1);
	BlockF = zeros(2*NU*Hu,NU*Hu);

	for i = 1:Hu
		for j = 1:NU
			% p.99のFをつくる． 1/下限上限 なので，0かどうかで場合分け
			if u_min(j) == 0
				F(2*NU*(i-1)+(2*j-1), NU*(i-1)+j) =  -1;
				f(2*NU*(i-1)+(2*j-1), 1)          =  0;
			else
				F(2*NU*(i-1)+(2*j-1), NU*(i-1)+j) =  -sign(u_min(j))/u_min(j);     % 下限
				f(2*NU*(i-1)+(2*j-1), 1)          =  sign(u_min(j));
			end
			if u_max(j) == 0
				F(2*NU*(i-1)+(2*j  ), NU*(i-1)+j) = 1;
				f(2*NU*(i-1)+(2*j  ), 1)          =  0;
			else
				F(2*NU*(i-1)+(2*j  ), NU*(i-1)+j) =   sign(u_max(j))/u_max(j);     % 上限
				f(2*NU*(i-1)+(2*j  ), 1)          = -sign(u_max(j));
			end
		end
	end
	for i = 1:Hp
		for j = 1:NX
			if x_min(j) == 0
				G(2*NX*(i-1)+(2*j-1), NX*(i-1)+j) =  -1;
				g(2*NX*(i-1)+(2*j-1), 1)          = 0;
			else
				G(2*NX*(i-1)+(2*j-1), NX*(i-1)+j) =  -sign(x_min(j))/x_min(j);     % 下限
				g(2*NX*(i-1)+(2*j-1), 1)          =  sign(x_min(j));
			end
			if x_max(j) == 0
				G(2*NX*(i-1)+(2*j  ), NX*(i-1)+j) = 1;
				g(2*NX*(i-1)+(2*j  ), 1)          =  0;
			else
				G(2*NX*(i-1)+(2*j  ), NX*(i-1)+j) =   sign(x_max(j))/x_max(j);     % 上限
				g(2*NX*(i-1)+(2*j  ), 1)          = -sign(x_max(j));
			end
		end
	end

	% BlockF
	temp = zeros(2*NU*Hu, NU);
	for i = 1:Hu
		temp = temp + F(:, NU*(i-1)+1:NU*i);
		BlockF(:, NU*(i-1)+1:NU*i) = temp;
	end

	flag = 1;
end

