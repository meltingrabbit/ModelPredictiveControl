function status = OutputFfsatPos_fromX(filename, x, ratio)
	% 2018-06-12
	% gnuplotで位置を描画するためのファイルを出力
	% filename			：出力ファイル名
	% x					：状態ベクトル
	% ratio				：倍率．省略可能


	if nargin == 2
		ratio = 1;
	end

	numMirSat = floor(length(x) / 6 - 1);

	status = -1;		% 異常終了


	str = '';

	% for i = 1:numMirSat
	% 	str = strcat(str, Vec2Str(ratio .* x(1:3)), '\r\n');
	% 	str = strcat(str, Vec2Str(ratio .* x(6*i+1:6*i+3)), '\r\n');
	% 	str = strcat(str, '\r\n\r\n');
	% end

	for i = 1:numMirSat
		str = strcat(str, Vec2Str(ratio .* x(6*i+1:6*i+3)), '\r\n');
	end
	str = strcat(str, Vec2Str(ratio .* x(6*1+1:6*1+3)), '\r\n');
	str = strcat(str, '\r\n\r\n');


	status = OutputFile(filename, str);

end


