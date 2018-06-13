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
	% cU(k)		: 0~Hu�܂ł�u�̗\���x�N�g��
	% ��U(k)	: 0~Hu�܂ł̇�u�̗\���x�N�g��

	% p.91
	% cZ(k) = Psi x(k) + Ups u(k-1) + Theta ��U(k)
	% cE(k) = cT(k) - Psi x(k) - Ups u(k-1)

	% ���� �� ��: ����ꎟ�s�����쐬
	% ######################################################


	%% �������Ȃ�
	persistent flag NX NU cQ cR Psi Ups Theta F f G g BlockF;
	% flag 		: ����Ăяo���p						: 1
	% NX 		: x�̃T�C�Y								: 1
	% NU 		: u�̃T�C�Y								: 1
	% c�͏��蕶�� mathcal
	% cQ 		: �d�ݍs��Q p.91 (3.3) 					: NX*Hp x NX*Hp
	% cR 		: �d�ݍs��R p.91 (3.4) 					: NU*Hu x NU*Hu
	% Psi 		: Psi p.91 (3.5), p.68					: NX*Hp x NX
	% Ups 		: Upsilon p.91 (3.5), p.68				: NX*Hp x NU
	% Theta 	: Theta p.91 (3.5), p.68				: NX*Hp x NU*Hu
	% F 	: cU�̐�������s�� p.99 (3.35) ��f������������ 		: 2*NU*Hu x NU*Hu ��������͏�������Ȃ̂�2*NU*Hu��
	% f 	: cU�̐�������x�N�g�� p.99 (3.35) ��f 				: 2*NU*Hu x 1
	% G 	: cZ�̐�������s�� p.99 (3.36) ��g������������ 		: 2*NX*Hp x NX*Hp ��������͏�������Ȃ̂�2*NX*Hp��
	% g 	: cZ�̐�������x�N�g�� p.99 (3.36) ��g 				: 2*NX*Hp x 1
	% BlockF: F�̗ݏ悵����� p.99�̉�����3�s�� 				: 2*NU*Hu x NU*Hu
	if isempty(flag)
		[flag, NX, NU, cQ, cR, Psi, Ups, Theta, F, f, G, g, BlockF] ...
					=  MPC_init(A, B, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max);
	end



	%% ����ꎟ�s��������
	% ���� �� ��: ����ꎟ�s�����쐬
	% Omega 	: �� p.100 (3.41)�̍��Ӎs��				:
	% omega 	: �� p.100 (3.41)�̉E�Ӄx�N�g��			:
	Omega = vertcat(BlockF, G*Theta);
	omega = vertcat(-BlockF(:,1:NU)*u0-f, -G*(Psi*x + Ups*u0)-g);


	%% �Q�ƋO������̎c�ԁC�Ǐ]�덷
	cT = repmat(r, Hp, 1);
	cE = cT - Psi*x - Ups*u0;					% �ÁDp.91 (3.6)

	cG = 2 * Theta.' * cQ * cE;					% G p.91 (3.11)
	cH = 2 * (Theta.' * cQ * Theta + cR);		% 2*H p.91 (3.12)

	% ������Ȃ��ƁCquadprog��
	% �x��: �w�b�V�A�����Ώ̂ł͂���܂���BH=(H+H')/2 �����Z�b�g���Ă��܂��B
	% ���ďo�Ă���D
	cH = round(cH, 10, 'significant');

	%% �œK���D�񎟌v��@
	opts = optimoptions(@quadprog,'Display','off');
	% opts = optimoptions(@quadprog,'Display','iter-detailed');
	% opts = optimoptions(@quadprog);
	dU = quadprog(cH, -cG.', Omega, omega, [], [], [], [], [], opts);
	du = dU(1:NU,1);

	u = u0 + du;
	% u = zeros(NU, 1);
end


%% �������֐�
% ���񂾂��Ă΂�邽�߁C�ω�����ϐ��͈����Ɏ��Ȃ�
function [flag, NX, NU, cQ, cR, Psi, Ups, Theta, F, f, G, g, BlockF] ...
		 = MPC_init(A, B, Hp, Hu, Q, Qf, R, x_min, x_max, u_min, u_max)
	% ����P�x�����Ă΂�鏀���֐�


	[NX, NU] = size(B);


	%% �d�ݍs�� p.91
	cQ = BlkdiagNTime(Q, Hp-1);
	cQ = blkdiag(cQ, Qf);
	cR = BlkdiagNTime(R, Hu);


	%% �\���s�� p.91
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
	for i = 1:Hu
		Theta( NX*(i-1)+1:NX*Hp , NU*(i-1)+1:NU*i ) = Ups( 1:NX*(Hp-i+1) , : );
	end


	%% ������� p.99
	F = zeros(2*NU*Hu,NU*Hu);
	f = zeros(2*NU*Hu,1);
	G = zeros(2*NX*Hp,NX*Hp);
	g = zeros(2*NX*Hp,1);
	BlockF = zeros(2*NU*Hu,NU*Hu);

	for i = 1:Hu
		for j = 1:NU
			% p.99��F������D 1/������� �Ȃ̂ŁC0���ǂ����ŏꍇ����
			if u_min(j) == 0
				F(2*NU*(i-1)+(2*j-1), NU*(i-1)+j) =  -1;
				f(2*NU*(i-1)+(2*j-1), 1)          =  0;
			else
				F(2*NU*(i-1)+(2*j-1), NU*(i-1)+j) =  -sign(u_min(j))/u_min(j);     % ����
				f(2*NU*(i-1)+(2*j-1), 1)          =  sign(u_min(j));
			end
			if u_max(j) == 0
				F(2*NU*(i-1)+(2*j  ), NU*(i-1)+j) = 1;
				f(2*NU*(i-1)+(2*j  ), 1)          =  0;
			else
				F(2*NU*(i-1)+(2*j  ), NU*(i-1)+j) =   sign(u_max(j))/u_max(j);     % ���
				f(2*NU*(i-1)+(2*j  ), 1)          = -sign(u_max(j));
			end
		end
	end
	for i = 1:Hp
		for j = 1:NX
			if x_min(j) == 0
				G(2*NX*(i-1)+(2*j-1), NX*(i-1)+j) =  -1;
				g(2*NX*(i-1)+(2*j-1), 1)          = 0;
			else
				G(2*NX*(i-1)+(2*j-1), NX*(i-1)+j) =  -sign(x_min(j))/x_min(j);     % ����
				g(2*NX*(i-1)+(2*j-1), 1)          =  sign(x_min(j));
			end
			if x_max(j) == 0
				G(2*NX*(i-1)+(2*j  ), NX*(i-1)+j) = 1;
				g(2*NX*(i-1)+(2*j  ), 1)          =  0;
			else
				G(2*NX*(i-1)+(2*j  ), NX*(i-1)+j) =   sign(x_max(j))/x_max(j);     % ���
				g(2*NX*(i-1)+(2*j  ), 1)          = -sign(x_max(j));
			end
		end
	end

	% BlockF
	temp = zeros(2*NU*Hu, NU);
	for i = 1:Hu
		temp = temp + F(:, NU*(i-1)+1:NU*i);
		BlockF(:, NU*(i-1)+1:NU*i) = temp;
	end

	flag = 1;
end

