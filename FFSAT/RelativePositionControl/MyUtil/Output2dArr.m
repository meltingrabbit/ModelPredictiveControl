function status = Output2dArr(filename, Arr, formatSpec)
	% 2018-06-13
	% ファイル出力
	% filename			：出力ファイル名
	% Arr				：配列
	% formatSpec		：フォーマット
	% https://jp.mathworks.com/help/matlab/ref/fprintf.html#btf98dm

	status = -1;		% 異常終了

	fileID = fopen(filename,'w');

	if fileID == -1
		return
	end

	fprintf(fileID, formatSpec ,Arr.');

	status = fclose(fileID);
end


