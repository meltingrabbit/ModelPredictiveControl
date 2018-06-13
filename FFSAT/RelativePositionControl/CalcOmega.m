function omega = CalcOmega(h, Constant)
	% 2018-06-12
	% ‹O“¹‚“x‚©‚ç‰~‹O“¹‚ÌŠp‘¬“x‚ğ‹‚ß‚éD

	omega = sqrt( Constant.Mu / (Constant.R + h)^3 );
end


