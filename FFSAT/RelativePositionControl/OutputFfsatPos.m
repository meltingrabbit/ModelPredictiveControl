function status = OutputFfsatPos(Ffsat, filename, ratio)
	% 2018-06-12
	% gnuplotで位置を描画するためのファイルを出力
	% Ffsat				：FFSAT構造体
	% filename			：出力ファイル名
	% ratio				：倍率．省略可能


	if nargin == 2
		ratio = 1;
	end


	status = -1;		% 異常終了


	str = '';

	for i = 1:Ffsat.Sat.MirSat.NUM
		str = strcat(str, Vec2Str(ratio .* Ffsat.Sat.ImgSat.x), '\r\n');
		str = strcat(str, Vec2Str(ratio .* Ffsat.Sat.MirSat.x{i}), '\r\n');
		str = strcat(str, '\r\n\r\n');
	end

	for i = 1:Ffsat.Sat.MirSat.NUM
		str = strcat(str, Vec2Str(ratio .* Ffsat.Sat.MirSat.x{i}), '\r\n');
	end
	str = strcat(str, Vec2Str(ratio .* Ffsat.Sat.MirSat.x{1}), '\r\n');
	str = strcat(str, '\r\n\r\n');


	status = OutputFile(filename, str);

end


