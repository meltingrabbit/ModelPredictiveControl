function [A,B] = GetHillEq(omega)
	% 2018-06-12
	% Hill方程式から状態空間モデルを生成する．
	% 座標系は軌道座標系
	% 軌道座標系は，衛星進行方向，軌道面外方向，地球中心方向をそれぞれ x,y,z軸とする直交座標系と定義される．

	A = zeros(6,6);
	B = zeros(6,3);

	I = eye(3);

	A(1:3, 4:6) = I;
	B(4:6, :)   = I;

	A(4,6) =  2*omega;
	A(6,4) = -2*omega;
	A(5,2) =  -omega^2;
	A(6,3) = 3*omega^2;


end