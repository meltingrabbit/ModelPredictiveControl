function omega = CalcOmega(h, Constant)
	% 2018-06-12
	% �O�����x����~�O���̊p���x�����߂�D

	omega = sqrt( Constant.Mu / (Constant.R + h)^3 );
end


