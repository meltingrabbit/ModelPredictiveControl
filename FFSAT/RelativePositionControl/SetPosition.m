function Ffsat = SetPosition(Ffsat, Constant, orientation)
	% 2018-06-12
	% 合成開口経，開口経，参照軌道から，各衛星の参照点からの相対位置を計算する．
	% Ffsat				：FFSAT構造体
	% orientation		：指向ベクトル（1x3）（大きさ不問）

	% 正規化
	orientation = orientation ./ norm(orientation);


	Ffsat.Sat.ImgSat.x = Ffsat.Sat.RefPoint.x + Ffsat.Telescope.F .* orientation;


	% orientationが軌道面垂直を向いていると死ぬw
	y_vec = Constant.ey;
	mirsat_plane_vec1 = cross(y_vec, orientation);
	mirsat_plane_vec2 = cross(orientation, mirsat_plane_vec1);


	for i = 1:Ffsat.Sat.MirSat.NUM
		rad = (i-1) * 2*pi / Ffsat.Sat.MirSat.NUM;
		d = ( Ffsat.Telescope.D - Ffsat.Sat.MirSat.D ) / 2;
		Ffsat.Sat.MirSat.x{i} = d .* ( mirsat_plane_vec1 .* cos(rad) + mirsat_plane_vec2 .* sin(rad) );
	end


end