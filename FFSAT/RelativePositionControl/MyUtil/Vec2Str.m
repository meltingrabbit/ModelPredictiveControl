function str = Vec2Str(x)
	% 2018-06-12
	% ベクトルをタブ区切り文字列に

	str = sprintf('%f',x(1));
	for i = 2:length(x)
		str = strcat(str, sprintf('\t%f',x(i)));
	end
end