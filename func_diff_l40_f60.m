function [ diff_mean40_mean60 ] = func_diff_l40_f60( delta )
%diff_l40_f60 Compute the difference of the mean of the last 40% elements
%with the mean of the first 60% elements of delta. This is to reduce the
%effect of the noise.
%   delta - only a portion of the full reference observation (full_delta)


length_first_60 = ceil(0.6*length(delta));
diff_mean40_mean60 = mean( delta((length_first_60+1):end) ) - mean( delta(1:length_first_60) );

end

