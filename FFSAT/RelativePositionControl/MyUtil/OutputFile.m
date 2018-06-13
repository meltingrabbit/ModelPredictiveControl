function status = OutputFile(filename, str, isApp)
	% 2018-06-12
	% �t�@�C���o��
	% filename			�F�o�̓t�@�C����
	% str				�F�o�͕�����
	% isApp				�F�ǋL���H�ȗ��\�D�f�t�H���g�l0


	if nargin == 2
		isApp = 0;
	end

	status = -1;		% �ُ�I��

	fileID = 0;			% �錾
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


