function status = Output2dArr(filename, Arr, formatSpec)
	% 2018-06-13
	% �t�@�C���o��
	% filename			�F�o�̓t�@�C����
	% Arr				�F�z��
	% formatSpec		�F�t�H�[�}�b�g
	% https://jp.mathworks.com/help/matlab/ref/fprintf.html#btf98dm

	status = -1;		% �ُ�I��

	fileID = fopen(filename,'w');

	if fileID == -1
		return
	end

	fprintf(fileID, formatSpec ,Arr.');

	status = fclose(fileID);
end


