function plan_fddt(input, t_samp, f_bin_width, min_drift, max_drift)

	nchan, ntime = size(input)
	# Block dimensions with units (Hz and seconds)
	t_block = ntime * t_samp
	f_block = nchan * f_bin_width
	# Min and max drifts across the block with units (Hz)
	min_freq_drift = ntime * t_samp * min_drift
	max_freq_drift = ntime * t_samp * max_drift

	# Minimum drift search resolution with integer frequency bin shifts (Hz/s)
	drift_resolution = f_block / t_block

	println("t_block: $t_block  f_block: $f_block  freq_drift: ($min_freq_drift, $max_freq_drift  drift_resolution: $drift_resolution")
	# TODO: Check maximum drift between two consecutive time samples to determine
	# initialization strategy (per FDMT paper at bottom of page 4)
	max_drift_per_samp = max_drift * t_samp
	println("Max freq drift per samp: $max_drift_per_samp
			Frequency bin width: $f_bin_width")

	# Convert delays to sample-space
	min_delay = Int(trunc(min_freq_drift / drift_resolution)) # TODO: Implement
	max_delay = Int(ceil( max_freq_drift / drift_resolution))

	nbatch = 1 # TODO: For GPU batching
	nsteps = Int(ceil(log2(ntime))) # Number of algorithm iterations

	# Create larger than necessary output size to fill to power of 2
	n_drift_trials = 2^nsteps
	out_size = (nchan, n_drift_trials)

	println("Out size: $out_size")
	output = zeros(Float32, out_size)

	println("Plan info:
			Input Size: ($nchan, $ntime)
			freq drift: ($min_freq_drift, $max_freq_drift)
			drift_resolution: $drift_resolution
			sample delays: ($min_delay, $max_delay)
			output_size: ($nchan, $n_drift_trials)")

	# TODO: Implement more complicated initialization with partial sums
    # Or just downsample along frequency axis
	if max_drift_per_samp > f_bin_width
		throw(error("Maximum frequency drift between two consecutive time samples is greater than the width of a frequency bin. Initialization of the output matrix requires additional work than the identity function."))
	# Else, the initial output data is simply the input data (Eq. 18 in the FDMT paper)
	else
		output[:,1:ntime] .= input
	end



	# Create the delays and source row data for later use in execution
	# TODO: See if it's better to do these calculations on the fly
	delays = Vector{Vector{Int32}}(undef, nsteps + 1) .= [[]]
	delays_mat = zeros(Int, (nsteps, n_drift_trials))
	# srcrows contains the indices of the rows used to create the current row of data
	srcrows = Vector{Vector{Tuple{Int32,Int32}}}(undef, nsteps + 1) .= [[]]
	srcrows_mat = zeros(Int, (nsteps, n_drift_trials,2))

	for step in 1:nsteps
		# Number of rows per subband of the current step
		nrow = 2^step
		resize!(srcrows[step], nrow)
		resize!(delays[step],  nrow)

		println("\nStep: $step  nrow: $nrow")

		stride = 2^(step - 1)

		for row in 1:nrow
			# srcrow0 = Int(floor(row / (2^(step-1)+1)) + 1)

			delays[step][row] = (row - 1) % 2

			if step == 1
				srcrows[step][row] = (1,2)
				delays[step][row] = row - 1
			else
				srcrow0 = Int(floor((row - 1) / 2)) + 1
				srcrows[step][row] = (srcrow0, srcrow0 + stride)
				# delays[step][row]  = row % stride == 0 ? 1 : 0 # TODO: Fix! This isn't right. Makes assumptions
				# delays[step][row]  = row % 2 == 0 ? 1 : 0 # TODO: Fix! This isn't right. Makes assumptions

			end

			# # Testing
			# delays_mat[step, row:nrow:end]     .=  delays[step][row]
			# srcrows_mat[step, row:nrow:end, 1] .=  srcrows[step][row][1]
			# srcrows_mat[step, row:nrow:end, 2] .=  srcrows[step][row][2]

			println("Srcrows: $(srcrows[step][row])\nDelays: $(delays[step][row])")
		end

	end

	return output, srcrows, delays
end

function fddt_exec(input, output, srcrows, delays)
	local nchan, ntime = size(input)
	local ndrift = size(output)[2]
	local nsteps = Int(log2(size(input)[2]))

	for step in 1:nsteps
		nrow = 2^step
		nsubband = Int(ntime / nrow)

		# Swap input and output
		if step > 1
			output, input = input, output
		end

		output = fddt_exec_step!(input,
						 output,
						 srcrows[step],
						 Int.(trunc.(collect(1:nrow) ./ 2)),
						 nsubband,
						 nrow,
						 step)
	end

	return output
end

function fddt_exec_step!(input, output, srcrows, delays, nsubband, nrow, step)
	nchan, ntime = size(input)
	# println("\nFDDT Exec Step: $step")
	# println("\tSrcrows: $srcrows  Delays: $delays")
	# println("\tNsubband: $nsubband  nrow: $nrow")

	for subband in 1:nsubband

		subband_offset = (subband - 1) * nrow
		# println("Subband offset: $subband_offset")

		for row in 1:nrow

			offset           = (row - 1) % 2^(step-1)
			delay            = delays[row]
			srcrow1, srcrow2 = srcrows[row]

			# println("Row Offset: $offset  Delay: $delay  Srcrows: ($srcrow1, $srcrow2)")

			for chan in 1:nchan

				outval = Float32(-1.0)

				if row - 1 + chan<= nchan
					# println("\tSB: $subband  row: $row  chan: $chan")
					# println("\tOut location: ($chan, $(subband_offset + row))")
					outval = input[chan, subband_offset + srcrow1] + input[chan + delay, subband_offset + srcrow2]
				end

				output[chan, subband_offset + row] = outval

			end
		end
	end

	# println("---FDDT Exec Step output: $output")
	return output
end

