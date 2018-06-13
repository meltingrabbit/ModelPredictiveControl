function B = BlkdiagNTime(A, n)
	% 2018-06-12
	% Blkdiag‚ğ“¯‚¶s—ñ‚Ån‰ñ‚â‚é
	% B = [A O O ....; O A O ....; .... ; .... O A ]


	ACell = repmat({A}, 1, n);
	B = blkdiag(ACell{:});
end


