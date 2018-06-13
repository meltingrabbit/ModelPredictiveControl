% FFSAT�̑��Έʒu��������f���\�������
% 2018-06-12

% close all
% clear
clear functions		% persistent�ϐ��̍폜

addpath('./MyUtil');
addpath('./Orbit');
addpath('./MPC');

%% �����ݒ�

Constant.R               = 6378100;			% �n�����a [m]
Constant.Mu              = 3.986004e14;		% �n�S�d�͒萔 [m^3/s^2]
Constant.P 	             = 4.617e-6;		% �n���ߖT�ł̑��z�t�ˈ� [N/m^2]
Constant.G               = 6.67408e-11;		% ���L���͒萔 [m^3/s^2/kg]
Constant.ex              = [1;0;0];			% �P�ʃx�N�g��
Constant.ey              = [0;1;0];			% �P�ʃx�N�g��
Constant.ez              = [0;0;1];			% �P�ʃx�N�g��

Orbit = GetOrbit(Constant);


% ������͑S��GEO���SHill���W�n
Ffsat.name               = 'FFSAT';
Ffsat.Telescope.F        = 10.0;			% �œ_���� [m]
Ffsat.Telescope.D        = 3.0;				% �����J���o [m]
Ffsat.Sat.RefPoint.x     = [0;0;0];			% FFSAT�̊���S�ʒu�i��ʒu�ix:�o�x����, y:�ܓx����, z:�n�S�����j����̃Y���j [m]
Ffsat.Sat.RefPoint.v     = [0;0;0];			% FFSAT�̊���S���x [m/s]
Ffsat.Sat.MirSat.NUM     = 6;				% ���q���� []
Ffsat.Sat.MirSat.MASS    = 50.0;			% ���q������ [kg]
Ffsat.Sat.MirSat.D       = 0.5;				% ���q���J���o [m]
Ffsat.Sat.MirSat.x       = cell(Ffsat.Sat.MirSat.NUM, 1);				% �ʒu�i����S�ʒu����̃Y���j [m]
Ffsat.Sat.MirSat.v       = cell(Ffsat.Sat.MirSat.NUM, 1);				% ���x�i����S�̑��Α��x�j [m/s]
for i = 1:Ffsat.Sat.MirSat.NUM
	Ffsat.Sat.MirSat.x{i} = [0;0;0];
	Ffsat.Sat.MirSat.v{i} = [0;0;0];
end
Ffsat.Sat.ImgSat.MASS    = 50.0;			% �B���q������ [kg]
Ffsat.Sat.ImgSat.x       = [0;0;0];			% �ʒu�i����S�ʒu����̃Y���j [m]
Ffsat.Sat.ImgSat.v       = [0;0;0];			% ���x�i����S�̑��Α��x�j [m/s]
% �X���X�^
Ffsat.Sat.MirSat.MAX_THRUST = 5.0e-6;				% �e���ő各�� [N]
Ffsat.Sat.ImgSat.MAX_THRUST = 5.0e-5;				% �e���ő各�� [N]
% �O��
Ffsat.Dist.C22           = @CalcDistC22;		% GEO�|�e���V�����o�x�����O��
Ffsat.Dist.Srp           = @CalcSrp;			% ���z�t�ˈ�
Ffsat.Dist.ThirdBody     = @CalcThirdBodyDist;	% ��O�V��



% Hill����
[HillA, HillB] = GetHillEq(Orbit.Geo.Omega);
% FFSAT�����ʒu
Ffsat = SetPosition(Ffsat, Constant, Constant.ez);


%% 7�@���̏�ԃx�N�g���C�s���`
x = vertcat(Ffsat.Sat.ImgSat.x, Ffsat.Sat.ImgSat.v);
for i = 1:Ffsat.Sat.MirSat.NUM
	x = vertcat(x, Ffsat.Sat.MirSat.x{i}, Ffsat.Sat.MirSat.v{i});
end
sizeX = 6*(Ffsat.Sat.MirSat.NUM+1);
sizeU = 3*(Ffsat.Sat.MirSat.NUM+1);
u = zeros(sizeU, 1);		% �������
d = zeros(sizeU, 1);		% �O��
A = BlkdiagNTime(HillA, Ffsat.Sat.MirSat.NUM+1);
B = BlkdiagNTime(HillB, Ffsat.Sat.MirSat.NUM+1);
% A = zeros(sizeX, sizeX);
% B = zeros(sizeX, sizeU);
% for i = 1:Ffsat.Sat.MirSat.NUM+1
% 	A( (i-1)*6+1:i*6, (i-1)*6+1:i*6 ) = HillA;
% 	B( (i-1)*6+1:i*6, (i-1)*3+1:i*3 ) = HillB;
% end



%% ������Ń����Q�N�b�^��
t0 = 0;
tn = 100000;
% tn = 1000;
dt = 10;
vt = t0:dt:tn;
count = 0;
countLog = 0;
% FORMAT_SPEC_LOG = strcat(repmat('%8e\t' , [1,sizeX+1]) , '\r\n');
FORMAT_SPEC_LOG = strcat(repmat('%f\t' , [1,sizeX+1]) , '\r\n');
% INTERVAL_LOG = 100;
INTERVAL_LOG = 10;
INTERVAL_MPC = 10;


% ��ԋ�ԃ��f��
sys = ss(A, B, zeros(1,sizeX), 0);
% ��ԋ�ԃ��f���i���U�j
sysd = c2d(sys, INTERVAL_MPC*dt);
Ad = sysd.A;
Bd = sysd.B;
% MPC�p�����^
r = x;			% �ڕW�l
Hp = 10;
Hu = 10;
Q = eye(sizeX);
Qf= Q;			% �Ō�����ɏd�݂�ς��Ȃ�
% R = zeros(sizeU,sizeU);
R = eye(sizeU);
x_min = repmat(-inf, [sizeX, 1]);
x_max = repmat( inf, [sizeX, 1]);
u_min = repmat(   0, [sizeX, 1]);
u_max = repmat( Ffsat.Sat.MirSat.MAX_THRUST / Ffsat.Sat.MirSat.MASS, [sizeX, 1]);
u_max(1:3) = repmat( Ffsat.Sat.ImgSat.MAX_THRUST / Ffsat.Sat.ImgSat.MASS, [3, 1]);


% ���ʊi�[�z��
output_tx = zeros(floor(length(vt) / INTERVAL_LOG)+1 , sizeX+1);
% �����Q�N�b�^
for t = vt
	% �O��
	d(1:3 ) = Ffsat.Dist.C22() + Ffsat.Dist.Srp(Orbit, Constant, t, 0.6, 0.25, Ffsat.Sat.ImgSat.MASS) + Ffsat.Dist.ThirdBody(Orbit, Constant, x(1:3), t);
	for i = 2:(Ffsat.Sat.MirSat.NUM+1)
		d( (i-1)*3+1:i*3 ) = Ffsat.Dist.C22() + Ffsat.Dist.Srp(Orbit, Constant, t, 0.6, 0.25, Ffsat.Sat.MirSat.MASS) + Ffsat.Dist.ThirdBody(Orbit, Constant, x((i-1)*6+1:i*6-3), t);

		% �m�F�p
		% aa = Ffsat.Dist.C22();
		% ab = Ffsat.Dist.Srp(Orbit, Constant, t, 0.6, 0.25, Ffsat.Sat.MirSat.MASS);
		% ac = Ffsat.Dist.ThirdBody(Orbit, Constant, x((i-1)*6+1:i*6-3), t);
	end

	% MPC
	if rem(count, INTERVAL_MPC) == 0
		u = MPC(x, u, Ad, Bd, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max);
	end

	% LOG
	if rem(count, INTERVAL_LOG) == 0
		countLog = countLog + 1;
		output_tx(countLog, :) = horzcat(t, x.');
		% OutputFfsatPos_fromX(strcat('./output/3Dpos/pos_', sprintf('%d', t) , '.dat'), x);
	end
	count = count + 1;

	dx = RungeKuttaSS(x, u+d, t, dt, A, B);
	x = x + dx;
end

Output2dArr('./output/tx.dat', output_tx, FORMAT_SPEC_LOG);

% OutputFfsatPos(Ffsat, 'FfsatPos.dat');





% ###############################
% �����Q�N�b�^�Ă���
% ###############################
% % x = vertcat(Ffsat.Sat.MirSat.x{2}, Ffsat.Sat.MirSat.v{2});
% x = vertcat(Ffsat.Sat.ImgSat.x, Ffsat.Sat.ImgSat.v);


% t0 = 0;
% tn = 100000;
% dt = 10;
% vt = t0:dt:tn;
% u = [0;0;0];
% dist = [0;0;0;0;0;0];

% OutputFile( 'x.dat', strcat(sprintf('%f',t0), '\t', Vec2Str(x), '\r\n' ), 0 );
% for t = vt
% 	dx = RungeKuttaSS(x, u, t, dt, HillA, HillB, dist);
% 	x = x + dx;
% 	if rem(t, 500) == 0
% 		OutputFile( 'x.dat', strcat(sprintf('%f',t), '\t', Vec2Str(x), '\r\n' ), 1 );
% 	end
% end
% ###############################







%% �O���v�Z�֐�
% ���ׂĉ����x��^����I�I
% a�����x�Cva�����x�x�N�g��
function va = CalcDistC22()
	% �o�x�����̂���

	a2dv = 365 * 24 * 60 * 60;		% �����x to ��V [m/s/yr]
	DeltaV = 1.314; 		% [m/s/yr]
	a = DeltaV / a2dv ;
	va = [-a;0;0];
end

function va = CalcSrp(Orbit, Constant, t, r, A, m)
	% ���z�t�ˈ����x�N�g���ŕԂ�
	% t 			: ����
	% r 			: ���ˌW��
	% A 			: ���˒f�ʐ�
	% m 			: �q���d��

	sunVec = Orbit.C.X2x(Orbit, t) * Orbit.VX.sun(Orbit, t);

	f    = Constant.P * r * A;
	va   = -(f/m) .* sunVec ./ norm(sunVec);
end

function va = CalcThirdBodyDist(Orbit, Constant, x, t)

	sunVecX  = Orbit.VX.sun(Orbit, t);
	moonVecX = Orbit.VX.moon(Orbit, t);
	% sunVec   = Orbit.C.X2x(Orbit, t) * sunVecX;
	% moonVec  = Orbit.C.X2x(Orbit, t) * moonVec;
	xVecX    = Orbit.VX.x0(Orbit, t) + Orbit.C.x2X(Orbit, t) * x;

	rSunVecX  = sunVecX - xVecX;
	rMoonVecX = moonVecX - xVecX;

	% vaMoon = Constant.G .* ( Orbit.Moon.MASS / norm(rMoonVecX)^3 .* rMoonVecX );
	vaMoon = [0;0;0];				% ���Ȃ�

	aSun   = Constant.G * Orbit.Sun.MASS / norm(rSunVecX)^2 - Orbit.Sun.CF;
	% aSun   = Constant.G * Orbit.Sun.MASS / norm(rSunVecX)^2 - norm(rSunVecX) * (2*pi/Orbit.Sun.T)^2;


	va = aSun .* rSunVecX ./ norm(rSunVecX) + vaMoon;
	va = Orbit.C.X2x(Orbit, t) * va;
end



