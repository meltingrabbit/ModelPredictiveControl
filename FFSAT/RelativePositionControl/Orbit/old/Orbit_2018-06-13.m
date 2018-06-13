function Orbit = GetOrbit(Orbit)
	% ��O�V��
	% GEO�z��

	%% �����ݒ�

	Orbit.name       = 'Orbit';

	Orbit.Geo.H      = 35786000;								% ��O�����xh_ref (GEO) [m]
	Orbit.Geo.Omega  = CalcOmega(Orbit.Geo.H, Constant);		% ����x�p���x (GEO) [rad/s]
	Orbit.Geo.T      = 2*pi / Orbit.Geo.Omega;					% ��O������ (GEO) [s]
	Orbit.Geo.R      = Orbit.Geo.H + Constant.R;				% ��O�����a (GEO) [m]

	Orbit.Sun.I      = deg2rad(23.44);				% �O���X�Ίp [rad]
	Orbit.Sun.MASS   = 1.9884e30;					% ���� [kg]
	Orbit.Sun.R      = 149597870700;				% �O�����a [m]
	Orbit.Moon.I     = deg2rad(5.14);				% �O���X�Ίp [rad]
	Orbit.Moon.MASS  = 7.347673e22;					% ���� [kg]
	Orbit.Moon.R     = 384400000;					% �O�����a [m]
	Orbit.Moon.T     = 29.5 * 24 * 3600;			% �i���z�n�����W�n�ł́j���� [s]

	% �ȉ��C�K���Ȑݒ�l
	% �ق�Ƃ͎��ςɂ���
	Orbit.Sun.Omega  = deg2rad(20);					% �n�����z�x�N�g�����瑪��������_�o�x [rad]
	Orbit.Moon.Omega = deg2rad(75);					% �n�����z�x�N�g�����瑪��������_�o�x [rad]



% GetSunVecX(Orbit, 100)
% GetSunVecX2(Orbit, 100)
% t = 60*60*3;
% xyz2XYZ(Orbit, [1;2;3], t)
% xyz2XYZ2(Orbit, [1;2;3], t)


end




function sunVecX = GetSunVecX(Orbit, t)
	% XYZ���W�n�ł̑��z�x�N�g��

	% �P�ʃx�N�g��
	% klm�͑��z�n�����W�n


	% XYZ���W�n�ł�klm���x�N�g��

	% klm���W�n�ł�XYZ���x�N�g��
	theta = pi/2 - Orbit.Sun.I;
	phi   = Orbit.Sun.Omega - pi/2;
	eX = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = pi/2;
	phi   = Orbit.Sun.Omega;
	eY = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = Orbit.Sun.I;
	phi   = Orbit.Sun.Omega + pi/2;
	eZ = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];


	C = vertcat(eX.', eY.', eZ.');			% ���W�ϊ��s��

	sunVeck = [1;0;0];						% klm�ł̑��z�x�N�g��

	sunVecX = C * sunVeck;
	sunVecX = sunVecX .* Orbit.Sun.R;
end

% % ���W�ϊ��s��g��Ȃ�Ver
% function sunVecX = GetSunVecX2(Orbit, t)
% 	% XYZ���W�n�ł̑��z�x�N�g��

% 	eY = [0;1;0]; 	% �P�ʃx�N�g��
% 	sunVecX = Ry(Orbit.Sun.I)*Rz(-Orbit.Sun.Omega)*eY;
% 	sunVecX = sunVecX .* Orbit.Sun.R;
% end



function moonVecX = GetMoonVecX(Orbit, t)
	% XYZ���W�n�ł̌��x�N�g��

	t0 = 0; 			% ������D�����ς���Ɗ�ʒu���ς��D
						% �܂�Ct0+t�ŉ��������(X����)�ɂ���̂ŁD

	% XYZ�͐Î~�O�����W�n
	% KLM�͌����W�n
	% klm�͑��z�n�����W�n


	% klm���W�n�ł�XYZ���x�N�g��






end



function X = xyz2XYZ(Orbit, x, t)
	% ���W�ϊ�
	% xzy��Hill���W�n
	% XYZ�͐Î~�O�����W�n
	% IMG_5448.JPG���m�[�g�Q��
	% �摜��XY���t

	t0 = 0; 			% ������D�����ς���Ɗ�ʒu���ς��D
						% �܂�Ct0+t�ŉ��������(X����)�ɂ���̂ŁD

	rad  = 2 * pi * (t0+t) / Orbit.Geo.T;

	% XYZ���W�n�ł�xyz���x�N�g��
	ex = [-sin(rad); cos(rad); 0];
	ey = [        0;        0;-1];
	ez = [-cos(rad);-sin(rad); 0];

	C = vertcat(ex.', ey.', ez.');			% ���W�ϊ��s��
	C = C.';								% XYZ�ōl���Ă���̂ŁC�]�u�i�t�s��j

	% �P�ʃx�N�g��
	eX = [1;0;0];
	eY = [0;1;0];
	% eZ = [0;0;1];

	X = C * x;
	% X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad)) + C * x;
end

% % ���W�ϊ��s��g��Ȃ�Ver
% function X = xyz2XYZ2(Orbit, x, t)
% 	% ���W�ϊ�
% 	% xzy��Hill���W�n
% 	% XYZ�͐Î~�O�����W�n
% 	% IMG_5448.JPG���m�[�g�Q��
% 	% �摜��XY���t

% 	t0 = 0; 			% ������D�����ς���Ɗ�ʒu���ς��D
% 						% �܂�Ct0+t�ŉ��������(X����)�ɂ���̂ŁD

% 	rad  = 2 * pi * (t0+t) / Orbit.Geo.T;

% 	% �P�ʃx�N�g��
% 	eX = [1;0;0];
% 	eY = [0;1;0];
% 	% eZ = [0;0;1];

% 	% X = Ry( -((pi/2)-rad)+pi )*Rx(pi/2)*x;
% 	% X = Rx(-pi/2)*Ry( ((pi/2)-rad)-pi )*x; 		% ��̋t�s��
% 	X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad)) + Rx(-pi/2)*Ry( ((pi/2)-rad)-pi )*x;
% end

