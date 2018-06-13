function status = OutputFfsatPos(Ffsat, filename, ratio)
	% 2018-06-12
	% gnuplot�ňʒu��`�悷�邽�߂̃t�@�C�����o��
	% Ffsat				�FFFSAT�\����
	% filename			�F�o�̓t�@�C����
	% ratio				�F�{���D�ȗ��\


	if nargin == 2
		ratio = 1;
	end


	status = -1;		% �ُ�I��


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


