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
Ffsat.Sat.MirSat.MAX_THRUST = 1.0e-3;				% �e���ő各�� [N]
Ffsat.Sat.ImgSat.MAX_THRUST = 1.0e-3;				% �e���ő各�� [N]
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


%% ���ΈʒuLOG�p
L_ImgMir = CalcDistance(Ffsat.Sat.MirSat.x{1}, Ffsat.Sat.ImgSat.x);				% �B���q��-���q���ԋ���
L_MirSid = CalcDistance(Ffsat.Sat.MirSat.x{1}, Ffsat.Sat.MirSat.x{2});			% �ׂ荇�������q���ԋ���
L_MirOpp = CalcDistance(Ffsat.Sat.MirSat.x{1}, Ffsat.Sat.MirSat.x{4});			% ���������������q���ԋ���


%% ����
SEED = 1;
rng(SEED);			% �����ݒ�


%% ������Ń����Q�N�b�^��
t0 = 0;
% tn = 10000;
% % tn = 1000;
% dt = 10;
tn = 1000;
dt = 1;
vt = t0:dt:tn;
count = 0;
countLog = 0;
% FORMAT_SPEC_LOG = strcat(repmat('%8e\t' , [1,sizeX+1]) , '\r\n');
% FORMAT_SPEC_LOG = strcat(repmat('%f\t' , [1,sizeX+1]) , '\r\n');
% INTERVAL_LOG = 100;
INTERVAL_LOG = 1;
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
R = zeros(sizeU,sizeU);
% R = eye(sizeU);
x_max = repmat( inf, [sizeX, 1]);
x_min = repmat(-inf, [sizeX, 1]);
u_max = repmat( Ffsat.Sat.MirSat.MAX_THRUST / Ffsat.Sat.MirSat.MASS, [sizeX, 1]);
u_max(1:3) = repmat( Ffsat.Sat.ImgSat.MAX_THRUST / Ffsat.Sat.ImgSat.MASS, [3, 1]);
% u_max = repmat( inf, [sizeX, 1]);
u_min = -u_max;


x_hat = x;		% �����MPC�́C�O��̗\�z���Ȃ��̂ō��̒l


% ���ʊi�[�z��
if INTERVAL_LOG == 1
	output_tx  = zeros(floor(length(vt) / INTERVAL_LOG) , sizeX+1);
	output_tu  = zeros(floor(length(vt) / INTERVAL_LOG) , sizeU+1);
	% output_err = zeros(floor(length(vt) / INTERVAL_LOG) , 8+1);
	output_err = zeros(floor(length(vt) / INTERVAL_LOG) , Ffsat.Sat.MirSat.NUM+1+1);
	output_errAll = zeros(floor(length(vt) / INTERVAL_LOG) , (Ffsat.Sat.MirSat.NUM+1)*4+1);
else
	output_tx  = zeros(floor(length(vt) / INTERVAL_LOG)+1 , sizeX+1);
	output_tu  = zeros(floor(length(vt) / INTERVAL_LOG)+1 , sizeU+1);
	% output_err = zeros(floor(length(vt) / INTERVAL_LOG)+1 , 8+1);
	output_err = zeros(floor(length(vt) / INTERVAL_LOG)+1 , Ffsat.Sat.MirSat.NUM+1+1);
	output_errAll = zeros(floor(length(vt) / INTERVAL_LOG)+1 , (Ffsat.Sat.MirSat.NUM+1)*4+1);
end
% �����Q�N�b�^
for t = vt
	% �O��
	d(1:3 ) = Ffsat.Dist.C22() + Ffsat.Dist.Srp(Orbit, Constant, t, 0.6, 0.25, Ffsat.Sat.ImgSat.MASS) + Ffsat.Dist.ThirdBody(Orbit, Constant, x(1:3), t);
	for i = 2:(Ffsat.Sat.MirSat.NUM+1)
		d( (i-1)*3+1:i*3 ) = Ffsat.Dist.C22() + Ffsat.Dist.Srp(Orbit, Constant, t, 0.6, 0.25, Ffsat.Sat.MirSat.MASS) + Ffsat.Dist.ThirdBody(Orbit, Constant, x((i-1)*6+1:i*6-3), t);

		% �m�F�p
		aa = Ffsat.Dist.C22();
		ab = Ffsat.Dist.Srp(Orbit, Constant, t, 0.6, 0.25, Ffsat.Sat.MirSat.MASS);
		ac = Ffsat.Dist.ThirdBody(Orbit, Constant, x((i-1)*6+1:i*6-3), t);
		na = norm(aa);
		nb = norm(ab);
		nc = norm(ac);
		aImg = HillA*x(1:6,:);
		aMir = HillA*x(13:18,:);
		nImg = norm(aImg(4:6,:));
		nMir = norm(aMir(4:6,:));
	end

	% MPC
	if rem(count, INTERVAL_MPC) == 0
		[u, x_hat] = MPC(x, x_hat, u, Ad, Bd, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max);
		t
	end

	% LOG
	if rem(count, INTERVAL_LOG) == 0
		countLog = countLog + 1;
		output_tx (countLog, :) = horzcat(t, x.');
		output_tu (countLog, :) = horzcat(t, u.');
		% output_err(countLog, :) = horzcat(t, GetRelPosErr(x, L_ImgMir, L_MirSid, L_MirOpp).');
		output_err(countLog, :) = horzcat(t, GetAbsPosErr(x, r).');
		output_errAll(countLog, :) = horzcat(t, GetAbsPosErrAll(x, r).');
		% OutputFfsatPos_fromX(strcat('./output/3Dpos/pos_', sprintf('%d', t) , '.dat'), x);
		% t
	end
	count = count + 1;

	dx = RungeKuttaSS(x, u+d, t, dt, A, B);
	x = x + dx;
end

[tmp, outputSize] = size(output_tx);
Output2dArr('./output/tx.dat', output_tx, strcat(repmat('%8e\t' , [1,outputSize]) , '\r\n'));
[tmp, outputSize] = size(output_tu);
Output2dArr('./output/tu.dat', output_tu, strcat(repmat('%8e\t' , [1,outputSize]) , '\r\n'));
[tmp, outputSize] = size(output_err);
Output2dArr('./output/err.dat', output_err, strcat(repmat('%8e\t' , [1,outputSize]) , '\r\n'));
[tmp, outputSize] = size(output_errAll);
Output2dArr('./output/errAll.dat', output_errAll, strcat(repmat('%8e\t' , [1,outputSize]) , '\r\n'));

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



%% LOG�̂��߂̊֐�
function e = GetRelPosErr(x, L_ImgMir, L_MirSid, L_MirOpp)
	% ���Έʒu�덷�z���Ԃ��D
	% Mir1-IMG, Mir1-Mir2, Mir1-Mir6, Mir1-Mir4, Mir3-IMG, Mir3-Mir4, Mir3-Mir2, Mir3-Mir6


	MIR_SAT_NUM = 6;
	e = zeros(4*2, 1);

	% IMG1, 3
	idx = 0;
	for i = [1,3]
		e(idx + 1) = CalcDistance(x(1:3), x(6*i+1:6*i+3)) - L_ImgMir;
		j = rem(i+1,MIR_SAT_NUM);
		if j==0
			j = MIR_SAT_NUM;
		end
		e(idx + 2) = CalcDistance(x(6*i+1:6*i+3), x(6*j+1:6*j+3)) - L_MirSid;
		j = rem(i-1,MIR_SAT_NUM);
		if j==0
			j = MIR_SAT_NUM;
		end
		e(idx + 3) = CalcDistance(x(6*i+1:6*i+3), x(6*j+1:6*j+3)) - L_MirSid;
		j = rem(i+3,MIR_SAT_NUM);
		if j==0
			j = MIR_SAT_NUM;
		end
		e(idx + 4) = CalcDistance(x(6*i+1:6*i+3), x(6*j+1:6*j+3)) - L_MirOpp;
		idx = idx + 4;
	end
end

function e = GetAbsPosErr(x, r)
	% ��Έʒu�덷�z���Ԃ��D

	MIR_SAT_NUM = 6;
	e = zeros(MIR_SAT_NUM+1, 1);

	for i = 1:(MIR_SAT_NUM+1)
		e(i) = CalcDistance(x(6*(i-1)+1:6*(i-1)+3), r(6*(i-1)+1:6*(i-1)+3));
	end
end

function e = GetAbsPosErrAll(x, r)
	% ��Έʒu�덷�z���Ԃ��D

	MIR_SAT_NUM = 6;
	e = zeros(MIR_SAT_NUM+1, 1);

	idx = 0;
	for i = 1:(MIR_SAT_NUM+1)
		e(idx+1) = CalcDistance(x(6*(i-1)+1:6*(i-1)+3), r(6*(i-1)+1:6*(i-1)+3));
		e(idx+2) = x(6*(i-1)+1) - r(6*(i-1)+1);
		e(idx+3) = x(6*(i-1)+2) - r(6*(i-1)+2);
		e(idx+4) = x(6*(i-1)+3) - r(6*(i-1)+3);
		idx = idx + 4;
	end
end
