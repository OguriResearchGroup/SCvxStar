%% interp_MEE: function description
function [x_his] = interp_MEE_from_cart(x0, xf, N, mu, options)
% *initial guess generation for direct trajectory optimization
% interpolate the initial and final state in MEE coordinates and 
% return an initial guess state history in the desired coordinates
	arguments
		x0 (6,1) double
		xf (6,1) double
		N (1,1) double
		mu (1,1) double
		options.ToF = []
		options.Nrev = []
	end
	if isempty(options.ToF) && isempty(options.Nrev)
		error('ToF or Nrev must be provided')
	end	

	ToF = options.ToF;
	Nrev = options.Nrev;

    import astro.conics.mee2cart
    import astro.conics.cart2mee

	mEOE0 = cart2mee(x0, mu);
	mEOEf = cart2mee(xf, mu);

	mEOE0(end) = wrapTo2Pi(mEOE0(end));
	while mEOEf(end) < mEOE0(end)
		mEOEf(end) = mEOEf(end) + 2*pi;
	end

	if isempty(Nrev)
		COE0vec = mEOE2COEvec(mEOE0);
		n0 		= sqrt(1/COE0vec(1)^3);
		Nrev 	= floor((n0 * ToF)/(2*pi));
	end

	mEOEf(end) = mEOEf(end) + Nrev * 2* pi;

	x_his = NaN(6,N);

	MEE_his = linspace_vec(mEOE0, mEOEf, N);
	for kk=1:N
		x_his(1:6,kk) = mee2cart(MEE_his(:,kk), mu);
	end

end
