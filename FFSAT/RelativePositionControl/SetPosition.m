function Ffsat = SetPosition(Ffsat, Constant, orientation)
	% 2018-06-12
	% �����J���o�C�J���o�C�Q�ƋO������C�e�q���̎Q�Ɠ_����̑��Έʒu���v�Z����D
	% Ffsat				�FFFSAT�\����
	% orientation		�F�w���x�N�g���i1x3�j�i�傫���s��j

	% ���K��
	orientation = orientation ./ norm(orientation);


	Ffsat.Sat.ImgSat.x = Ffsat.Sat.RefPoint.x + Ffsat.Telescope.F .* orientation;


	% orientation���O���ʐ����������Ă���Ǝ���w
	y_vec = Constant.ey;
	mirsat_plane_vec1 = cross(y_vec, orientation);
	mirsat_plane_vec2 = cross(orientation, mirsat_plane_vec1);


	for i = 1:Ffsat.Sat.MirSat.NUM
		rad = (i-1) * 2*pi / Ffsat.Sat.MirSat.NUM;
		d = ( Ffsat.Telescope.D - Ffsat.Sat.MirSat.D ) / 2;
		Ffsat.Sat.MirSat.x{i} = d .* ( mirsat_plane_vec1 .* cos(rad) + mirsat_plane_vec2 .* sin(rad) );
	end


end