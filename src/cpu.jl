export fddt


"""
	function fddt(input::AbstractArray, t_samp::AbstractFloat, f_bin_width::AbstractFloat,
				  min_drift::AbstractFloat, max_drift::AbstractFloat; use_gpu=false)

Dedoppler the input data using a Taylor-Tree-based algorithm

Inputs:
		input:       Input data array to dedoppler
		t_samp:      Time sample duration in seconds
		f_bin_width: Frequency channel width in Hertz
		min_drift:   Minimum drift rate to search in Hertz per second
		max_drift:   Maximum drift rate to search in Hertz per second
		use_gpu:     Whether to dispatch to a CUDA-enabled GPU or not

Outputs the dedopplered array
"""
function fddt(input::AbstractArray, t_samp::AbstractFloat, f_bin_width::AbstractFloat,
			  min_drift::AbstractFloat, max_drift::AbstractFloat; use_gpu=false)

	buffer_1, buffer_2, freq_scrunch_ratio = create_plan(input, t_samp, f_bin_width,
														 min_drift, max_drift)
	# If dispatching to the GPU, check for functionality and push data to GPU
	if use_gpu
		@assert CUDA.functional(true)
		input = cu(input)
		buffer_1 = cu(buffer_1)
		buffer_2 = cu(buffer_2)
	end

	fddt_out = exec(input, buffer_1, buffer_2, freq_scrunch_ratio)

	# Return GPU data back to the host before returning
	if use_gpu
		fddt_out = Array(fddt_out)
	end


	# TODO: Output custom struct with metadata
	return fddt_out
end

"""
	function create_plan(input, t_samp, f_bin_width, min_drift, max_drift)

Setup and create one-time information for FDDT execution
"""
function create_plan(input, t_samp, f_bin_width, min_drift, max_drift)

	# TODO: Valid input parameter checks

	nchan, ntime = size(input)
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

"""
	function init!(input::AbstractArray, buffer_1::AbstractArray,
						freq_scrunch_ratio::Integer, input_ntime::Integer, nchan::Integer)

Initialize FDDT data for execution

Reduce (sum) along frequency dimension by freq_scrunch_ratio. Initialize only the ntime
first drift rate rows with the input data. This effectively creates zero padding below
ntime to create power of 2 length.
"""
function init!(input::AbstractArray, buffer_1::AbstractArray,
					freq_scrunch_ratio::Integer, input_ntime::Integer, nchan::Integer)
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
	return nothing
end

function exec(input::AbstractArray, buffer_1::AbstractArray, buffer_2::AbstractArray,
			  freq_scrunch_ratio::Integer)
	input_ntime = size(input)[2]
	nchan, ndrift = size(buffer_1)
	nsteps = Int(log2(ndrift))
	srcrows = Vector{Tuple{Int32, Int32}}(undef, ndrift)

	init!(input, buffer_1, freq_scrunch_ratio, input_ntime, nchan)

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

		exec_step!(buffer_1,
				   buffer_2,
				   srcrows,
				   delays,
				   nsubband,
				   nrow,
				   nchan)
	end

	return buffer_2
end

"""
	function fddt_exec_step!(input::AbstractArray, output::AbstractArray, srcrows, delays,
							 nsubband, nrow, nchan)

Execute a single step of the FDDT algorithm
"""
function exec_step!(input::AbstractArray, output::AbstractArray, srcrows, delays, nsubband, nrow, nchan)

	for subband in 1:nsubband

		subband_offset = (subband - 1) * nrow

		for row in 1:nrow

			delay            = delays[row]
			srcrow1, srcrow2 = srcrows[row]

			for chan in 1:nchan

				outval = Float32(-1.0)

				if row - 1 + chan <= nchan
					outval = input[chan,         subband_offset + srcrow1] +
							 input[chan + delay, subband_offset + srcrow2]
				end

				output[chan, subband_offset + row] = outval

			end
		end
	end
	return nothing
end
