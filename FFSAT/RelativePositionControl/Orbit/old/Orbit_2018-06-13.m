function Orbit = GetOrbit(Orbit)
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
	Orbit.Moon.I     = deg2rad(5.14);				% 軌道傾斜角 [rad]
	Orbit.Moon.MASS  = 7.347673e22;					% 質量 [kg]
	Orbit.Moon.R     = 384400000;					% 軌道半径 [m]
	Orbit.Moon.T     = 29.5 * 24 * 3600;			% （太陽地球座標系での）周期 [s]

	% 以下，適当な設定値
	% ほんとは時変にする
	Orbit.Sun.Omega  = deg2rad(20);					% 地球太陽ベクトルから測った昇交点経度 [rad]
	Orbit.Moon.Omega = deg2rad(75);					% 地球太陽ベクトルから測った昇交点経度 [rad]



% GetSunVecX(Orbit, 100)
% GetSunVecX2(Orbit, 100)
% t = 60*60*3;
% xyz2XYZ(Orbit, [1;2;3], t)
% xyz2XYZ2(Orbit, [1;2;3], t)


end




function sunVecX = GetSunVecX(Orbit, t)
	% XYZ座標系での太陽ベクトル

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

	sunVeck = [1;0;0];						% klmでの太陽ベクトル

	sunVecX = C * sunVeck;
	sunVecX = sunVecX .* Orbit.Sun.R;
end

% % 座標変換行列使わないVer
% function sunVecX = GetSunVecX2(Orbit, t)
% 	% XYZ座標系での太陽ベクトル

% 	eY = [0;1;0]; 	% 単位ベクトル
% 	sunVecX = Ry(Orbit.Sun.I)*Rz(-Orbit.Sun.Omega)*eY;
% 	sunVecX = sunVecX .* Orbit.Sun.R;
% end



function moonVecX = GetMoonVecX(Orbit, t)
	% XYZ座標系での月ベクトル

	t0 = 0; 			% 基準時刻．これを変えると基準位置が変わる．
						% つまり，t0+tで黄道から上(X軸上)にいるので．

	% XYZは静止軌道座標系
	% KLMは月座標系
	% klmは太陽地球座標系


	% klm座標系でのXYZ基底ベクトル






end



function X = xyz2XYZ(Orbit, x, t)
	% 座標変換
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

	% 単位ベクトル
	eX = [1;0;0];
	eY = [0;1;0];
	% eZ = [0;0;1];

	X = C * x;
	% X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad)) + C * x;
end

% % 座標変換行列使わないVer
% function X = xyz2XYZ2(Orbit, x, t)
% 	% 座標変換
% 	% xzyはHill座標系
% 	% XYZは静止軌道座標系
% 	% IMG_5448.JPGかノート参照
% 	% 画像はXYが逆

% 	t0 = 0; 			% 基準時刻．これを変えると基準位置が変わる．
% 						% つまり，t0+tで黄道から上(X軸上)にいるので．

% 	rad  = 2 * pi * (t0+t) / Orbit.Geo.T;

% 	% 単位ベクトル
% 	eX = [1;0;0];
% 	eY = [0;1;0];
% 	% eZ = [0;0;1];

% 	% X = Ry( -((pi/2)-rad)+pi )*Rx(pi/2)*x;
% 	% X = Rx(-pi/2)*Ry( ((pi/2)-rad)-pi )*x; 		% 上の逆行列
% 	X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad)) + Rx(-pi/2)*Ry( ((pi/2)-rad)-pi )*x;
% end

