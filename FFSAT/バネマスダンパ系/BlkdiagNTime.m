function B = BlkdiagNTime(A, n)
	% 2018-06-12
	% Blkdiag�𓯂��s���n����
	% B = [A O O ....; O A O ....; .... ; .... O A ]


	ACell = repmat({A}, 1, n);
	B = blkdiag(ACell{:});
end


