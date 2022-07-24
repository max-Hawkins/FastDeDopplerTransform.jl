### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ b50f3c28-09ed-11ed-08f2-37407421f5b5
begin
	using Pkg
	Pkg.activate("..")

	using JLD
	using Revise
	using Plots
	using BenchmarkTools

	using FastDeDopplerTransform
end

# ╔═╡ 525f3f0f-6d4e-4925-996c-342ab3c943db
pwd()

# ╔═╡ 2899618d-0b5f-47dc-9653-49756762e5a6
begin
	_file = load("../voyager_test.jld")
	small_data = _file["raw_data"]
	fch1     = _file["fch1"]
	t_samp   = _file["tsamp"]
	f_bin_width       = _file["df"] * -1e6

	nchan, ntime = size(small_data)
	ndrift = ntime
	
	nothing
end

# ╔═╡ 94f90bfd-f4bf-45b6-b5e3-b0e6fe385c25
f_bin_width

# ╔═╡ 338660f8-c953-400c-aaab-4d263e01647d
t_samp

# ╔═╡ 6d9d2516-4498-4619-b712-6ea29d7b80cd
f_bin_width / t_samp

# ╔═╡ 1ee14061-6236-4456-bc5f-9a3c7ee2927c
_t = Int(ceil(4.0 / (f_bin_width / t_samp)))

# ╔═╡ 639289f7-d1f6-4b2b-94ec-00fc3326a854
begin
	# Find closest (but greater) integer to nchan that's divisible by freq_scrunch_ratio
	
end

# ╔═╡ 4f687826-b8df-4ee8-984d-14bdab92fda7
2^20 / (_t + 5)

# ╔═╡ 5f7cc3f1-94c7-4f76-aafd-f84f62b7ab93
2^20 % (_t+5)

# ╔═╡ 85a94335-77ca-415f-996b-b80fe41dc5c4
0.152 * ntime * t_samp / f_bin_width

# ╔═╡ 7bb11fd0-71b5-40f0-affd-116f1b228e1a
begin
	min_drift = 0.0   # Hz/s
	max_drift = 0.5 # Hz/s
	
	buffer_1, buffer_2, srcrows, freq_scrunch_ratio = FastDeDopplerTransform.plan_fddt(small_data, t_samp, f_bin_width, min_drift, max_drift)

	fddt_small_data = FastDeDopplerTransform.fddt_exec(small_data, buffer_1, buffer_2, srcrows, freq_scrunch_ratio)

	heatmap(fddt_small_data', 
		title = "Voyager 1 FDDT Output - Drifts: ($min_drift, $max_drift)",
        xlabel = "Frequency channel",
        ylabel = "Drift Rate (Normalized)",
		yflip=true)
end

# ╔═╡ 64fe80f7-45bb-4f07-9e63-7a6a99ee6ef5
begin
	in = zeros(Int32, (5,4))
	for i in 1:size(in)[2]
		in[:,i] .= collect(1:size(in)[1])
	end
end

# ╔═╡ 9eb9a01d-173c-4c71-ba08-083d99bb3871
heatmap(in',yflip=true)

# ╔═╡ e1f19d81-d747-4254-903f-f95adf3e0897
begin
	scrunch = 2
	# out_nchan = Int(floor(size(in)[1] / scrunch))
	# reshaped = reshape(in[1:out_nchan*scrunch,:], (scrunch, out_nchan  ,size(in)[2]))
end

# ╔═╡ 1501b9f5-f62b-4665-8abb-7f700bf851d0
size(in)

# ╔═╡ 878afce6-6bef-4726-9a9e-f5d9943bd99e
sum(reshaped, dims=2)

# ╔═╡ 009e81d3-cfd6-4aff-9a86-5a76c30b3d2b
begin
	in_nchan, in_ntime = size(in)
	out_nchan, out_trials = Int(floor(in_nchan / scrunch)), in_ntime

	_out = zeros(Int, (out_nchan, out_trials))
	for t in 1:in_ntime
		for c in 1:out_nchan
			in_start_chan = Int((c-1) * scrunch + 1)
			println(typeof(in_start_chan))
			_out[c, t] = sum(in[in_start_chan:in_start_chan + scrunch - 1, t])
		end
	end
	_out
end

# ╔═╡ 2cc80dd9-03e7-4f9a-a9d0-091489b20a28
heatmap(_out', yflip=true)

# ╔═╡ Cell order:
# ╠═b50f3c28-09ed-11ed-08f2-37407421f5b5
# ╠═525f3f0f-6d4e-4925-996c-342ab3c943db
# ╠═2899618d-0b5f-47dc-9653-49756762e5a6
# ╠═94f90bfd-f4bf-45b6-b5e3-b0e6fe385c25
# ╠═338660f8-c953-400c-aaab-4d263e01647d
# ╠═6d9d2516-4498-4619-b712-6ea29d7b80cd
# ╠═1ee14061-6236-4456-bc5f-9a3c7ee2927c
# ╠═639289f7-d1f6-4b2b-94ec-00fc3326a854
# ╠═4f687826-b8df-4ee8-984d-14bdab92fda7
# ╠═5f7cc3f1-94c7-4f76-aafd-f84f62b7ab93
# ╠═85a94335-77ca-415f-996b-b80fe41dc5c4
# ╠═7bb11fd0-71b5-40f0-affd-116f1b228e1a
# ╠═64fe80f7-45bb-4f07-9e63-7a6a99ee6ef5
# ╠═9eb9a01d-173c-4c71-ba08-083d99bb3871
# ╠═e1f19d81-d747-4254-903f-f95adf3e0897
# ╠═1501b9f5-f62b-4665-8abb-7f700bf851d0
# ╠═878afce6-6bef-4726-9a9e-f5d9943bd99e
# ╠═009e81d3-cfd6-4aff-9a86-5a76c30b3d2b
# ╠═2cc80dd9-03e7-4f9a-a9d0-091489b20a28
