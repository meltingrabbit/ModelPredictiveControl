function u = MPC(x, u0, A, B, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max)
	% ���f���\������
	% u 		: ���̃X�e�b�v�̏o��					: NU
	% x 		: ��ԕϐ�								: NX
	% u0 		: 1�X�e�b�v�O�̓���						: NU
	% A 		: ���U�� ��ԕ�����						: NX x NX
	% B 		: ���U�� ��ԕ�����						: NX x NU
	% r 		: �ڕW��ԁi���͌Œ��z�肵�Ă���j	: NX
	% Hp		: �\���z���C�Y��						: 1
	% Hu		: ����z���C�Y��						: 1
	% Q			: ��ԗʂɑ΂���d�ݍs�� 				: NX x NX
	% Qf		: ��ԗʂɑ΂���ŏI�d�ݍs��			: NX x NX
	% R			: ���͕ω��ɑ΂���d�ݍs��				: NU x NU
	% x_min		: �������								: NX
	% x_max		: �������								: NX
	% u_min		: �������								: NU
	% u_max		: �������								: NU


	% ######################################################
	% MPC
	% ######################################################
	% c�͏��蕶�� mathcal
	% x(k+1) = A x(k) + B u(k)
	% y(k+1) = Cy x(k)				: ����x�N�g��
	% z(k+1) = Cz x(k)				: �o�̓x�N�g��
	% r(k) 							: �ڕW�l�C�Q�ƋO��

	% p.90
	% cZ(k) 	: 1~Hp�܂ł�z�i�����ł� =y =x�j�̗\���x�N�g��
	% cT(k) 	: 1~Hp�܂ł�r�̖ڕW�l�x�N�g��
	% ��U(k)	: 0~Hu�܂ł̇�u�̗\���x�N�g��

	% p.91
	% cZ(k) = Pis x(k) + Ups u(k-1) + Theta ��U(k)
	% cE(k) = cT(k) - Pis x(k) - Ups u(k-1)
	% ######################################################


	%% �������Ȃ�
	persistent flag NX NU cQ cR;
	% flag 		: ����Ăяo���p						: 1
	% NX 		: x�̃T�C�Y								: 1
	% NU 		: u�̃T�C�Y								: 1
	% c�͏��蕶�� mathcal
	% cQ 		: �d�ݍs��Q p.91 (3.3) 					: NX*Hp x NX*Hp
	% cR 		: �d�ݍs��R p.91 (3.4) 					: NU*Hu x NU*Hu
	% Psi 		: Psi p.91 (3.5), p68					: NX*Hp x NX
	% Ups 		: Upsilon p.91 (3.5), p68				: NX*Hp x NU
	% Theta 	: Theta p.91 (3.5), p68					: NX*Hp x NU*Hu

	if isempty(flag)
		[flag, NX, NU, cQ, cR] ...
					=  MPC_init(x, u0, A, B, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max);
	end



	u = zeros(NU, 1);


end


function [flag, NX, NU, cQ, cR] = MPC_init(x, u0, A, B, r, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max)
	% ����P�x�����Ă΂�鏀���֐�


	[NX, NU] = size(B);
	cQ = BlkdiagNTime(Q, Hp-1);
	cQ = blkdiag(cQ, Qf);
	cR = BlkdiagNTime(R, Hu);

	Psi   = zeros(NX*Hp, NX);
	Ups   = zeros(NX*Hp, NU);
	Theta = zeros(NX*Hp, NU*Hu);

	temp_Apow    = eye(NX);
	temp_UpsPre  = zeros(NX, NU);
	for i = 1:Hp
		Ups(NX*(i-1)+1:NX*i, :) = temp_Apow * B + temp_UpsPre;
		temp_UpsPre = Ups(NX*(i-1)+1:NX*i, :);

		temp_Apow = temp_Apow * A;
		Psi(NX*(i-1)+1:NX*i, :) = temp_Apow;
	end



	flag = 1;
end

