function status = OutputFfsatPos_fromX(filename, x, ratio)
	% 2018-06-12
	% gnuplot�ňʒu��`�悷�邽�߂̃t�@�C�����o��
	% filename			�F�o�̓t�@�C����
	% x					�F��ԃx�N�g��
	% ratio				�F�{���D�ȗ��\


	if nargin == 2
		ratio = 1;
	end

	numMirSat = floor(length(x) / 6 - 1);

	status = -1;		% �ُ�I��


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


