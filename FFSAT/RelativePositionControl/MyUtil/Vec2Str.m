function str = Vec2Str(x)
	% 2018-06-12
	% �x�N�g�����^�u��؂蕶�����

	str = sprintf('%f',x(1));
	for i = 2:length(x)
		str = strcat(str, sprintf('\t%f',x(i)));
	end
end