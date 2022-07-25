function plan_fddt(input, t_samp, f_bin_width, min_drift, max_drift)

	nchan, ntime = size(input)
	# Block dimensions with units
	t_block = ntime * t_samp # seconds
	f_block = nchan * f_bin_width # Hz
	# Min and max drifts across the block with units
	min_freq_drift = ntime * t_samp * min_drift # Hz
	max_freq_drift = ntime * t_samp * max_drift # Hz

	# Minimum drift search resolution with integer frequency bin shifts (Hz/s)
	# Also doubles as the greatest drift rate possible without frequency reduction
	drift_resolution = f_bin_width / t_samp # f_block / t_block

	# Frequency reduction ratio
	freq_scrunch_ratio = Int(ceil(max_drift / drift_resolution))

	# Update drift_resolution to account for frequency scrunching
	# this will increase f_bin_width
	drift_resolution *= freq_scrunch_ratio

	# Update number of output frequency channels to account for rounding
	# Will throw away at most (freq_scrunch_ratio - 1) channels. TODO: Is this acceptable?
	out_nchan = Int(floor(nchan / freq_scrunch_ratio))

	# Convert delays to sample-space (accounting for frequency scrunching)
	min_delay = Int(trunc(min_freq_drift / (f_bin_width * freq_scrunch_ratio))) # TODO: Implement
	max_delay = Int(ceil( max_freq_drift / (f_bin_width * freq_scrunch_ratio)))

	# Check that max delay isn't greater than the number of frequency channels
	if max_delay >= out_nchan
		throw(error("Maximum number of channel shifts to meet requested drift rate \
		 rate exceeds number of frequency channels. Reduce max drift rate.
		 Requested max drift rate: $max_drift
		 Sample Delays (with frequency scrunching factor of $freq_scrunch_ratio): $max_delay
		 Output number of channels: $out_nchan"))
	end

	# Create larger than necessary output size to fill to power of 2
	# Reduce along frequency dimension if the max drift rate is too high (freq_scrunch_ratio)
	nsteps = Int(ceil(log2(max_delay))) # Number of algorithm iterations
	n_drift_trials = 2^nsteps
	out_size = (out_nchan, n_drift_trials)

	# Allocate and return two working buffers to not overwrite input data
	buffer_1 = zeros(Float32, out_size)
	buffer_2 = zeros(Float32, out_size)

	return buffer_1, buffer_2, freq_scrunch_ratio
end

function fddt_exec(input, buffer_1, buffer_2, freq_scrunch_ratio)
	input_ntime = size(input)[2]
	nchan, ndrift = size(buffer_1)
	nsteps = Int(log2(size(buffer_1)[2]))
	srcrows = Vector{Tuple{Int32, Int32}}(undef, ndrift)

	# Initialize data
		# Reduce (sum) along frequency dimension by freq_scrunch_ratio
			# TODO: Do with reshaping and views/slices?
		# Initialize only the ntime first drift rate rows with the input data
			# This effectively creates zero padding below ntime to create power of 2 length
	if freq_scrunch_ratio > 1
		for t in 1:input_ntime
			for c in 1:nchan-freq_scrunch_ratio
				in_start_chan = Int((c-1) * freq_scrunch_ratio + 1)
				buffer_1[c, t] = sum(input[in_start_chan:in_start_chan + freq_scrunch_ratio - 1, t])
			end
		end
	else
		buffer_1[:,1:input_ntime] .= input
	end

	for step in 1:nsteps
		nrow = 2^step
		stride = 2^(step-1)
		nsubband = Int(ndrift / nrow)

		# Calculate delays
		delays = Int.(trunc.(collect(1:nrow) ./ 2))

		# Calculate srcrows
		for r in 1:stride
			for i in 0:1
				srcrows[(r-1)*2+i+1] = (r, r+stride)
			end
		end

		# Swap buffers
		if step > 1
			buffer_2, buffer_1 = buffer_1, buffer_2
		end

		buffer_2 = fddt_exec_step!(buffer_1,
						 buffer_2,
						 srcrows,
						 delays,
						 nsubband,
						 nrow,
						 nchan)
	end

	return buffer_2
end

function fddt_exec_step!(input, output, srcrows, delays, nsubband, nrow, nchan)

	for subband in 1:nsubband

		subband_offset = (subband - 1) * nrow

		for row in 1:nrow

			delay            = delays[row]
			srcrow1, srcrow2 = srcrows[row]

			for chan in 1:nchan

				outval = Float32(-1.0)

				if row - 1 + chan <= nchan
					outval = input[chan, subband_offset + srcrow1] + input[chan + delay, subband_offset + srcrow2]
				end

				output[chan, subband_offset + row] = outval

			end
		end
	end
	return output
end

