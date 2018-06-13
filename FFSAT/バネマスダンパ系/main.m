% �o�l�}�X�_���p�n��MPC�e�X�g
% 2018-06-12

% close all
clear
clear functions		% persistent�ϐ��̍폜


%% physical parameter
m = 1.0;
c = 1.0;
k = 0.8;


sizeX = 2;
sizeU = 1;
A = [0, 1; -k/m, -c/m];
B = [0; 1/m];

x = [30;0];
u = 0;		% �������
d = 0;		% �O��




t0 = 0;
tn = 10;
% tn = 1000;
dt = 0.1;
vt = t0:dt:tn;
count = 0;
countLog = 0;
FORMAT_SPEC_LOG = strcat(repmat('%f\t' , [1,sizeX+1]) , '\r\n');
INTERVAL_LOG = 1;
INTERVAL_MPC = 10;


% ��ԋ�ԃ��f��
sys = ss(A, B, zeros(1,sizeX), 0);
% ��ԋ�ԃ��f���i���U�j
sysd = c2d(sys, INTERVAL_MPC*dt);
Ad = sysd.A;
Bd = sysd.B;
% MPC�p�����^
r = [-5;0];			% �ڕW�l
Hp = 10;
Hu = 5;
Q = eye(sizeX);
Qf= Q;			% �Ō�����ɏd�݂�ς��Ȃ�
% R = zeros(sizeU,sizeU);
R = eye(sizeU);
x_min = [-inf; -inf];
x_max = repmat( inf, [sizeX, 1]);
u_min = [-2];
u_max = repmat( inf, [sizeX, 1]);


% ���ʊi�[�z��
if INTERVAL_LOG == 1
	output_tx = zeros(floor(length(vt) / INTERVAL_LOG) , sizeX+1);
	output_tu = zeros(floor(length(vt) / INTERVAL_LOG) , sizeU+1);
else
	output_tx = zeros(floor(length(vt) / INTERVAL_LOG)+1 , sizeX+1);
	output_tu = zeros(floor(length(vt) / INTERVAL_LOG)+1 , sizeU+1);
end
% �����Q�N�b�^
for t = vt
	% �O��
	% �Ȃ�

	% MPC
	if rem(count, INTERVAL_MPC) == 0
		u = MPC(x, u, Ad, Bd, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max);
	end

	% LOG
	if rem(count, INTERVAL_LOG) == 0
		countLog = countLog + 1;
		output_tx(countLog, :) = horzcat(t, x.');
		output_tu(countLog, :) = horzcat(t, u.');
	end
	count = count + 1;

	dx = RungeKuttaSS(x, u+d, t, dt, A, B);
	x = x + dx;
end


figure(1);
subplot(2,1,1);
plot(output_tx(:,1), output_tx(:,2));
subplot(2,1,2);
plot(output_tx(:,1), output_tx(:,3));

figure(2);
plot(output_tu(:,1), output_tu(:,2));


% Output2dArr('./output/tx.dat', output_tx, FORMAT_SPEC_LOG);






