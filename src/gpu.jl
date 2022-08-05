"""
	function exec_step_kernel(input, output, nsubband, nrow, nchan, step)

Execute a single step of the FDDT algorithm dispatched to a CUDA-enabled GPU
"""
function exec_step_kernel(input, output, nsubband, nrow, nchan, step)

	# Calculate indices
	chan = threadIdx().x + (blockIdx().z - 1) * blockDim().x
	row  = blockIdx().y

	nrow = nrow[1]
	step = step[1]
	nsubband = nsubband[1]
	nchan = nchan[1]

	# Calculate delay
	delay = Int(trunc(((row-1) % (2^step) + 1) / 2))

	# Exit early if invalid thread
	if chan + delay > nchan
		return nothing
	end

	# Calculate source rows to sum together
	stride  = 2 ^ (step - 1)
	subband = floor((row - 1) / (nrow / nsubband))
	srcrow1 = Int(ceil(row / 2) + subband * stride)
	srcrow2 = srcrow1 + stride

	# Write to output
	output[chan, row] = input[chan,         srcrow1] +
						input[chan + delay, srcrow2]

	return nothing
end

"""
	function exec(input::CuArray, buffer_1::CuArray, buffer_2::CuArray, freq_scrunch_ratio::Integer)

Execute the entire FDDT algorithm on a CUDA-enabled GPU
"""
function exec(input::CuArray, buffer_1::CuArray, buffer_2::CuArray, freq_scrunch_ratio::Integer)
	input_ntime = size(input)[2]
	nchan, ndrift = size(buffer_1)
	nsteps = Int(log2(ndrift))

	# Initialize data
	# Reduce (sum) along frequency dimension by freq_scrunch_ratio
	# Initialize only the ntime first drift rate rows with the input data
		# This effectively creates zero padding below ntime to create power of 2 length
	if freq_scrunch_ratio == 1
		# If no scrunching, fill buffer directly with input data
		CUDA.@sync buffer_1[:,1:input_ntime] .= input
	else
		# Else sum every freq_scrunch_ratio channels together
		# This was faster than a naive kernel. TODO: Evaluated tuned kernel performance
		input_reshaped = reshape(input[begin:nchan*freq_scrunch_ratio,:], freq_scrunch_ratio, nchan, input_ntime)
		buffer_1_sum_dest = reshape(@view(buffer_1[:,begin:input_ntime]), 1, nchan, input_ntime);
		CUDA.@sync sum!(buffer_1_sum_dest, input_reshaped)
	end

	# Calculate launch configuration
	_threads = (min(1024, nchan))
	_blocks  = (1, ndrift, Int(ceil(nchan / _threads[1])))

	# For each step in the FDDT algorithm, launch the execution kernel
	for step in 1:nsteps
		nrow = 2^step
		nsubband = Int(ndrift / nrow)

		# Swap buffers
		if step > 1
			buffer_2, buffer_1 = buffer_1, buffer_2
		end

		d_nsubband = cu([nsubband])
		d_ndrift   = cu([ndrift])
		d_nchan    = cu([nchan])
		d_step     = cu([step])

		# Execute kernel per step
		CUDA.@sync @cuda threads=_threads blocks=_blocks exec_step_kernel(buffer_1,
						 buffer_2,
						 d_nsubband,
						 d_ndrift,
						 d_nchan,
						 d_step)
	end

	# Return the output - always buffer_2
	return buffer_2
end