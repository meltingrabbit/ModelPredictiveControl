function status = OutputFile(filename, str, isApp)
	% 2018-06-12
	% ファイル出力
	% filename			：出力ファイル名
	% str				：出力文字列
	% isApp				：追記か？省略可能．デフォルト値0


	if nargin == 2
		isApp = 0;
	end

	status = -1;		% 異常終了

	fileID = 0;			% 宣言
	if isApp
		fileID = fopen(filename,'a');
	else
		fileID = fopen(filename,'w');
	end

	if fileID == -1
		return
	end

	fprintf(fileID, str);
	status = fclose(fileID);
end


