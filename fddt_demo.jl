### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ b50f3c28-09ed-11ed-08f2-37407421f5b5
begin
	using Pkg
	Pkg.activate("..")

	using Downloads
	using HDF5
	using Revise
	using Plots
	using CUDA
	using BenchmarkTools

	using FastDeDopplerTransform
end

# ╔═╡ 525f3f0f-6d4e-4925-996c-342ab3c943db
begin
	voyager_url = "http://blpd14.ssl.berkeley.edu/voyager_2020/single_coarse_channel/single_coarse_guppi_59046_80036_DIAG_VOYAGER-1_0011.rawspec.0000.h5"
	h5 = h5open(Downloads.download(voyager_url))
	full_spectrogram = h5["data"][:,1,:]
	freq_range = range(659935, length=150)
	spectrogram = full_spectrogram[freq_range,:]

	fch1        = h5["data"]["fch1"][]
	t_samp      = h5["data"]["tsamp"][]
	f_bin_width = h5["data"]["foff"][] * -1e6

	min_drift = 0.0 # Hz/s
	max_drift = 0.5 # Hz/s
	
	size(spectrogram), size(full_spectrogram)
end

# ╔═╡ 7bb11fd0-71b5-40f0-affd-116f1b228e1a
begin
	local buffer_1, buffer_2, srcrows, freq_scrunch_ratio = FastDeDopplerTransform.plan_fddt(spectrogram, t_samp, f_bin_width, min_drift, max_drift)

	fddt_small_data = FastDeDopplerTransform.fddt_exec(spectrogram, buffer_1, buffer_2, srcrows, freq_scrunch_ratio)
	
	nothing
end

# ╔═╡ 3673612f-e0e4-44f2-9e0b-89dbddd52ca1
heatmap(fddt_small_data', 
		title = "Voyager 1 FDDT Output - Drifts: ($min_drift, $max_drift)",
        xlabel = "Frequency channel",
        ylabel = "Drift Rate (Normalized)",
		yflip=true)

# ╔═╡ 64fe80f7-45bb-4f07-9e63-7a6a99ee6ef5
begin
	local buffer_1, buffer_2, freq_scrunch_ratio = FastDeDopplerTransform.plan_fddt(full_spectrogram, t_samp, f_bin_width, min_drift, max_drift)

	fddt_data = FastDeDopplerTransform.fddt_exec(full_spectrogram, buffer_1, buffer_2, freq_scrunch_ratio)

	nothing
end

# ╔═╡ 04aff514-1942-4547-8b34-1fe5809127a3
heatmap(fddt_data', 
		title = "Voyager 1 FDDT Output - Drifts: ($min_drift, $max_drift)",
        xlabel = "Frequency channel",
        ylabel = "Drift Rate (Normalized)",
		yflip=true)

# ╔═╡ 2d3169cb-3c78-42a6-9367-0536466b8478


# ╔═╡ 57939149-de57-4c05-bece-c5b3f720e90a
begin
	max_idx = findmax(fddt_data)[2][1]
	heatmap(fddt_data[max_idx-100:max_idx+100,:]', 
		title = "Voyager 1 FDDT Output - Drifts: ($min_drift, $max_drift)",
        xlabel = "Frequency channel",
        ylabel = "Drift Rate (Normalized)",
		yflip=true)

end

# ╔═╡ 36cd79b4-bd90-423e-902c-411a8b2ecad2
md"""
## View Top Hits
 TODO: Slide to pick N-th top hit
"""

# ╔═╡ 7e3b3eee-96c0-4889-92cc-0d21db6542e3
begin
	local max_hit = 15

	local temp_data = similar(fddt_data)
	temp_data .= fddt_data
	local max_idx = findmax(temp_data)[2][1]
	println("Max idx: $(findmax(temp_data))")

	# Zero out max_hit - 1 top hits
	for hit in 1:max_hit-1
		temp_data[max_idx-25:max_idx+25,:] .= Float32(0.0)
		
		max_idx = findmax(temp_data)[2][1]
		println("Max idx: $(findmax(temp_data))")
	end

	
	heatmap(temp_data[max_idx-100:max_idx+100,:]', 
		title = "Voyager 1 FDDT Output - Drifts: ($min_drift, $max_drift)",
        xlabel = "Frequency channel",
        ylabel = "Drift Rate (Normalized)",
		yflip=true)

end

# ╔═╡ 8f9a8174-bc78-4e1b-9a19-9b4d89647e9f
md"""
# Benchmarking
"""

# ╔═╡ e01dc75c-b155-4ddc-8f83-a0dd55f77b78
md"""
## CPU - Small Spectrogram
"""

# ╔═╡ c4aa5c98-63c5-4d23-b16a-148698db5836
begin
	local buffer_1, buffer_2, freq_scrunch_ratio = FastDeDopplerTransform.plan_fddt(spectrogram, t_samp, f_bin_width, min_drift, max_drift)

	 @benchmark FastDeDopplerTransform.fddt_exec($spectrogram, $buffer_1, $buffer_2, $freq_scrunch_ratio)
end

# ╔═╡ 85072289-1ca0-44f2-b84f-5dda3faebb02
md"""
## CPU - Full Spectrogram
"""

# ╔═╡ 1cffcac3-f59e-4e3d-856e-fea1beedc3d9
begin
	local buffer_1, buffer_2, freq_scrunch_ratio = FastDeDopplerTransform.plan_fddt(full_spectrogram, t_samp, f_bin_width, min_drift, max_drift)

	@benchmark FastDeDopplerTransform.fddt_exec($full_spectrogram, $buffer_1, $buffer_2, $freq_scrunch_ratio)
end

# ╔═╡ ee571ece-2947-4675-a3be-019776155127
md"""
## GPU - Small Spectrogram
"""

# ╔═╡ 48ac2ae8-1773-4873-a131-ae7cfa47a667
begin
	local buffer_1, buffer_2, freq_scrunch_ratio = FastDeDopplerTransform.plan_fddt(spectrogram, t_samp, f_bin_width, min_drift, max_drift)

	local d_spectrogram = cu(spectrogram)
	local d_buffer_1 = cu(buffer_1)
	local d_buffer_2 = cu(buffer_2)
	local d_freq_scrunch_ratio = cu(freq_scrunch_ratio)

	CUDA.allowscalar(true)

	@benchmark CUDA.@sync FastDeDopplerTransform.fddt_exec($d_spectrogram, $d_buffer_1, $d_buffer_2, $d_freq_scrunch_ratio)
end

# ╔═╡ fc8248fe-482f-4f87-915d-c6ff396166ea
md"""
## GPU - Full Spectrogram
"""

# ╔═╡ 0327f759-0375-45d0-b9d7-b72b7ac7fd1d
begin
	local buffer_1, buffer_2, freq_scrunch_ratio = FastDeDopplerTransform.plan_fddt(full_spectrogram, t_samp, f_bin_width, min_drift, max_drift)

	local d_full_spectrogram = cu(full_spectrogram)
	local d_buffer_1 = cu(buffer_1)
	local d_buffer_2 = cu(buffer_2)
	local d_freq_scrunch_ratio = cu(freq_scrunch_ratio)

	CUDA.allowscalar(true)

	CUDA.@sync FastDeDopplerTransform.fddt_exec(d_full_spectrogram, d_buffer_1, d_buffer_2, d_freq_scrunch_ratio)
end

# ╔═╡ Cell order:
# ╠═b50f3c28-09ed-11ed-08f2-37407421f5b5
# ╠═525f3f0f-6d4e-4925-996c-342ab3c943db
# ╠═7bb11fd0-71b5-40f0-affd-116f1b228e1a
# ╠═3673612f-e0e4-44f2-9e0b-89dbddd52ca1
# ╠═64fe80f7-45bb-4f07-9e63-7a6a99ee6ef5
# ╠═04aff514-1942-4547-8b34-1fe5809127a3
# ╠═2d3169cb-3c78-42a6-9367-0536466b8478
# ╠═57939149-de57-4c05-bece-c5b3f720e90a
# ╟─36cd79b4-bd90-423e-902c-411a8b2ecad2
# ╠═7e3b3eee-96c0-4889-92cc-0d21db6542e3
# ╟─8f9a8174-bc78-4e1b-9a19-9b4d89647e9f
# ╟─e01dc75c-b155-4ddc-8f83-a0dd55f77b78
# ╠═c4aa5c98-63c5-4d23-b16a-148698db5836
# ╟─85072289-1ca0-44f2-b84f-5dda3faebb02
# ╠═1cffcac3-f59e-4e3d-856e-fea1beedc3d9
# ╟─ee571ece-2947-4675-a3be-019776155127
# ╠═48ac2ae8-1773-4873-a131-ae7cfa47a667
# ╟─fc8248fe-482f-4f87-915d-c6ff396166ea
# ╠═0327f759-0375-45d0-b9d7-b72b7ac7fd1d
