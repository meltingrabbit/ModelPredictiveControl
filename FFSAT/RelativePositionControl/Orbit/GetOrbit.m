function Orbit = GetOrbit(Constant)
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
	Orbit.Sun.T      = 365.256363004 * 24 * 3600;	% ���� [s]
	Orbit.Sun.CF     = Orbit.Sun.R * (2*pi/Orbit.Sun.T)^2;		% ���z�n�����W�n�ł̉��S�́i�����x�j centrifugal force [m/s^2]
	Orbit.Moon.I     = deg2rad(5.14);				% �O���X�Ίp [rad]
	Orbit.Moon.MASS  = 7.347673e22;					% ���� [kg]
	Orbit.Moon.R     = 384400000;					% �O�����a [m]
	Orbit.Moon.T     = 29.5 * 24 * 3600;			% �i���z�n�����W�n�ł́j���� [s]

	% �ȉ��C�K���Ȑݒ�l
	% �ق�Ƃ͎��ςɂ���
	Orbit.Sun.Omega  = deg2rad(20);					% �n�����z�x�N�g�����瑪��������_�o�x [rad]
	Orbit.Moon.Omega = deg2rad(75);					% �n�����z�x�N�g�����瑪��������_�o�x [rad]

	Orbit.C.x2X      = @x2X;						% Hill���W�n to GEO���ʍ��W�n x2X(Orbit, t)
	Orbit.C.X2x      = @X2x;
	Orbit.C.k2X      = @k2X;						% �n�����z���W�ϊ� to GEO���ʍ��W�n k2X(Orbit, t)
	Orbit.C.X2k      = @X2k;
	Orbit.C.k2K      = @k2K;						% �n�����z���W�ϊ� to ���O���ʍ��W�n k2K(Orbit, t)
	Orbit.C.K2k      = @K2k;
	Orbit.VX.x0      = @x0onX;						% GEO���ʍ��W�n����݂�Hill���W�n���_�x�N�g��
	Orbit.VX.sun     = @GetSunVecX;					% GEO���ʍ��W�n����݂����z�x�N�g��
	Orbit.VX.moon    = @GetMoonVecX;				% GEO���ʍ��W�n����݂����x�N�g��


% GetSunVecX(Orbit, 100)
% GetSunVecX2(Orbit, 100)
% t = 60*60*3;
% xyz2XYZ(Orbit, [1;2;3], t)
% xyz2XYZ2(Orbit, [1;2;3], t)

% GetMoonVecX(Orbit, 0)
% GetMoonVecX(Orbit, 24*60*60*10)


end




%% �x�N�g���擾�֐�
function X = x0onX(Orbit, t)
	% GEO���ʍ��W�n����݂�Hill���W�n���_�x�N�g��

	rad  = 2 * pi * t / Orbit.Geo.T;

	% �P�ʃx�N�g��
	eX = [1;0;0];
	eY = [0;1;0];
	% eZ = [0;0;1];

	X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad));
end

function sunVecX = GetSunVecX(Orbit, t)
	% XYZ���W�n�ł̑��z�x�N�g��

	sunVeck = [1;0;0];						% klm�ł̑��z�x�N�g��

	sunVecX = Orbit.C.k2X(Orbit, t) * sunVeck;
	sunVecX = sunVecX .* Orbit.Sun.R;
end

function moonVecX = GetMoonVecX(Orbit, t)
	% XYZ���W�n�ł̌��x�N�g��

	t0 = 0; 			% ������D�����ς���Ɗ�ʒu���ς��D
						% �܂�Ct0+t�ŉ��������(X����)�ɂ���̂ŁD

	% XYZ�͐Î~�O�����W�n
	% KLM�͌����W�n
	% klm�͑��z�n�����W�n

	% �P�ʃx�N�g��
	eK = [1;0;0];
	eL = [0;1;0];
	% eM = [0;0;1];

	rad  = 2 * pi * (t0+t) / Orbit.Moon.T;
	moonVecK = Orbit.Moon.R .* (eK .* cos(rad) + eL .* sin(rad));						% KLM�ł̌��x�N�g��

	moonVecX = Orbit.C.k2X(Orbit, t) * Orbit.C.K2k(Orbit, t) * moonVecK;
end


% % Hill��x����X��
% function X = xyz2XYZ(Orbit, x, t)
% 	% �P�ʃx�N�g��
% 	eX = [1;0;0];
% 	eY = [0;1;0];
% 	% eZ = [0;0;1];

% 	X = Orbit.C.x2X(Orbit, t) * x;
% 	% X = Orbit.Geo.R .* (eX .* cos(rad) + eY .* sin(rad)) + Orbit.C.x2X(Orbit, t) * x;
% end



%% ���W�ϊ�

function C = x2X(Orbit, t)
	% ���W�ϊ�
	% Hill���W�n to GEO���ʍ��W�n x2X(Orbit, t)

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
end
function C = X2x(Orbit, t)
	C = x2X(Orbit, t).';
end

function C = k2X(Orbit, t)
	% �n�����z���W�ϊ� to GEO���ʍ��W�n k2X(Orbit, t)

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
end
function C = X2k(Orbit, t)
	C = k2X(Orbit, t).';
end

function C = k2K(Orbit, t)
	% �n�����z���W�ϊ� to ���O���ʍ��W�n

	% klm���W�n�ł�KLM���x�N�g��
	theta = pi/2 - Orbit.Moon.I;
	phi   = Orbit.Moon.Omega - pi/2;
	eK = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = pi/2;
	phi   = Orbit.Moon.Omega;
	eL = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];
	theta = Orbit.Moon.I;
	phi   = Orbit.Moon.Omega + pi/2;
	eM = [sin(theta)*cos(phi);sin(theta)*sin(phi);cos(theta)];

	C = vertcat(eK.', eL.', eM.');			% ���W�ϊ��s��
end
function C = K2k(Orbit, t)
	C = k2K(Orbit, t).';
end

