function [A,B] = GetHillEq(omega)
	% 2018-06-12
	% Hill�����������ԋ�ԃ��f���𐶐�����D
	% ���W�n�͋O�����W�n
	% �O�����W�n�́C�q���i�s�����C�O���ʊO�����C�n�����S���������ꂼ�� x,y,z���Ƃ��钼�����W�n�ƒ�`�����D

	A = zeros(6,6);
	B = zeros(6,3);

	I = eye(3);

	A(1:3, 4:6) = I;
	B(4:6, :)   = I;

	A(4,6) =  2*omega;
	A(6,4) = -2*omega;
	A(5,2) =  -omega^2;
	A(6,3) = 3*omega^2;


end