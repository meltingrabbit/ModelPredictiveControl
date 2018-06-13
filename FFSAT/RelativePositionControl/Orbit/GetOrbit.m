function Orbit = GetOrbit(Constant)
	% 第三天体
	% GEO想定

	%% 初期設定

	Orbit.name       = 'Orbit';

	Orbit.Geo.H      = 35786000;								% 基準軌道高度h_ref (GEO) [m]
	Orbit.Geo.Omega  = CalcOmega(Orbit.Geo.H, Constant);		% 基準高度角速度 (GEO) [rad/s]
	Orbit.Geo.T      = 2*pi / Orbit.Geo.Omega;					% 基準軌道周期 (GEO) [s]
	Orbit.Geo.R      = Orbit.Geo.H + Constant.R;				% 基準軌道半径 (GEO) [m]

	Orbit.Sun.I      = deg2rad(23.44);				% 軌道傾斜角 [rad]
	Orbit.Sun.MASS   = 1.9884e30;					% 質量 [kg]
	Orbit.Sun.R      = 149597870700;				% 軌道半径 [m]
	Orbit.Sun.T      = 365.256363004 * 24 * 3600;	% 周期 [s]
	Orbit.Sun.CF     = Orbit.Sun.R * (2*pi/Orbit.Sun.T)^2;		% 太陽地球座標系での遠心力（加速度） centrifugal force [m/s^2]
	Orbit.Moon.I     = deg2rad(5.14);				% 軌道傾斜角 [rad]
	Orbit.Moon.MASS  = 7.347673e22;					% 質量 [kg]
	Orbit.Moon.R     = 384400000;					% 軌道半径 [m]
	Orbit.Moon.T     = 29.5 * 24 * 3600;			% （太陽地球座標系での）周期 [s]

	% 以下，適当な設定値
	% ほんとは時変にする
	Orbit.Sun.Omega  = deg2rad(20);					% 地球太陽ベクトルから測った昇交点経度 [rad]
	Orbit.Moon.Omega = deg2rad(75);					% 地球太陽ベクトルから測った昇交点経度 [rad]

	Orbit.C.x2X      = @x2X;						% Hill座標系 to GEO平面座標系 x2X(Orbit, t)
	Orbit.C.X2x      = @X2x;
	Orbit.C.k2X      = @k2X;						% 地球太陽座標変換 to GEO平面座標系 k2X(Orbit, t)
	Orbit.C.X2k      = @X2k;
	Orbit.C.k2K      = @k2K;						% 地球太陽座標変換 to 月軌道面座標系 k2K(Orbit, t)
	Orbit.C.K2k      = @K2k;
	Orbit.VX.x0      = @x0onX;						% GEO平面座標系からみたHill座標系原点ベクトル
	Orbit.VX.sun     = @GetSunVecX;					% GEO平面座標系からみた太陽ベクトル
	Orbit.VX.moon    = @GetMoonVecX;				% GEO平面座標系からみた月ベクトル


% GetSunVecX(Orbit, 100)
% GetSunVecX2(Orbit, 100)
% t = 60*60*3;
% xyz2XYZ(Orbit, [1;2;3], t)
% xyz2XYZ2(Orbit, [1;2;3], t)

% GetMoonVecX(Orbit, 0)
% GetMoonVecX(Orbit, 24*60*60*10)


end




%% ベクトル取得関数
function X = x0onX(Orbit, t)
	% GEO平面座標系からみたHill座標系原点ベクトル

	rad  = 2 * pi * t / Orbit.Geo.T;

	% 単位ベクトル
	eX = [1;0;0];
	eY = [0;1;0];
	% eZ = [0;0;1];

	X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad));
end

function sunVecX = GetSunVecX(Orbit, t)
	% XYZ座標系での太陽ベクトル

	sunVeck = [1;0;0];						% klmでの太陽ベクトル

	sunVecX = Orbit.C.k2X(Orbit, t) * sunVeck;
	sunVecX = sunVecX .* Orbit.Sun.R;
end

function moonVecX = GetMoonVecX(Orbit, t)
	% XYZ座標系での月ベクトル

	t0 = 0; 			% 基準時刻．これを変えると基準位置が変わる．
						% つまり，t0+tで黄道から上(X軸上)にいるので．

	% XYZは静止軌道座標系
	% KLMは月座標系
	% klmは太陽地球座標系

	% 単位ベクトル
	eK = [1;0;0];
	eL = [0;1;0];
	% eM = [0;0;1];

	rad  = 2 * pi * (t0+t) / Orbit.Moon.T;
	moonVecK = Orbit.Moon.R .* (eK .* cos(rad) + eL .* sin(rad));						% KLMでの月ベクトル

	moonVecX = Orbit.C.k2X(Orbit, t) * Orbit.C.K2k(Orbit, t) * moonVecK;
end


% % HillのxからXに
% function X = xyz2XYZ(Orbit, x, t)
% 	% 単位ベクトル
% 	eX = [1;0;0];
% 	eY = [0;1;0];
% 	% eZ = [0;0;1];

% 	X = Orbit.C.x2X(Orbit, t) * x;
% 	% X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad)) + Orbit.C.x2X(Orbit, t) * x;
% end



%% 座標変換

function C = x2X(Orbit, t)
	% 座標変換
	% Hill座標系 to GEO平面座標系 x2X(Orbit, t)

	% xzyはHill座標系
	% XYZは静止軌道座標系
	% IMG_5448.JPGかノート参照
	% 画像はXYが逆

	t0 = 0; 			% 基準時刻．これを変えると基準位置が変わる．
						% つまり，t0+tで黄道から上(X軸上)にいるので．

	rad  = 2 * pi * (t0+t) / Orbit.Geo.T;

	% XYZ座標系でのxyz基底ベクトル
	ex = [-sin(rad); cos(rad); 0];
	ey = [        0;        0;-1];
	ez = [-cos(rad);-sin(rad); 0];

	C = vertcat(ex.', ey.', ez.');			% 座標変換行列
	C = C.';								% XYZで考えているので，転置（逆行列）
end
function C = X2x(Orbit, t)
	C = x2X(Orbit, t).';
end

function C = k2X(Orbit, t)
	% 地球太陽座標変換 to GEO平面座標系 k2X(Orbit, t)

	% 単位ベクトル
	% klmは太陽地球座標系

	% XYZ座標系でのklm基底ベクトル

	% klm座標系でのXYZ基底ベクトル
	theta = pi/2 - Orbit.Sun.I;
	phi   = Orbit.Sun.Omega - pi/2;
	eX = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = pi/2;
	phi   = Orbit.Sun.Omega;
	eY = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = Orbit.Sun.I;
	phi   = Orbit.Sun.Omega + pi/2;
	eZ = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];

	C = vertcat(eX.', eY.', eZ.');			% 座標変換行列
end
function C = X2k(Orbit, t)
	C = k2X(Orbit, t).';
end

function C = k2K(Orbit, t)
	% 地球太陽座標変換 to 月軌道面座標系

	% klm座標系でのKLM基底ベクトル
	theta = pi/2 - Orbit.Moon.I;
	phi   = Orbit.Moon.Omega - pi/2;
	eK = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = pi/2;
	phi   = Orbit.Moon.Omega;
	eL = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = Orbit.Moon.I;
	phi   = Orbit.Moon.Omega + pi/2;
	eM = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];

	C = vertcat(eK.', eL.', eM.');			% 座標変換行列
end
function C = K2k(Orbit, t)
	C = k2K(Orbit, t).';
end

