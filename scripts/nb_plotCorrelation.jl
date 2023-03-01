### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 6086b32e-a3ce-11ed-336c-1dd69a0fed42
begin
	using Pkg
	Pkg.activate("../")
	Pkg.add(["JSON3","CairoMakie","AlgebraOfGraphics","DataFrames","StatsBase","XLSX","RCall","MAT","ColorSchemes","CategoricalArrays","ColorTypes","DataFramesMeta"])
	Pkg.develop("MakieThemes")
	using JSON3
	#using MakieCore
	#using Makie
	using MakieThemes

	using CairoMakie
	using AlgebraOfGraphics
	using DataFrames
	using StatsBase
	using XLSX
	using RCall
	using MAT
	using ColorSchemes
	using CategoricalArrays
	using DataFramesMeta
using ColorTypes
	
end

# ╔═╡ af862c8f-9d52-4d29-85fa-4289dae0c6d5
CairoMakie.set_theme!(MakieThemes.ggthemr(:fresh))

# ╔═╡ 4fde6495-cb2a-4855-af11-36ac865cd0d9
md"""
## Loading Decoding
"""

# ╔═╡ 19eb5619-c2d1-47d1-9c09-2484c9e9f00b
function read_json!(dList,t,s,crossOrWithin)

		for split = (crossOrWithin == "Cross" ? 0 : 1:10)
			if split == 0
				
			jsonpath = "SubjectAnalysis_$crossOrWithin/deep/$t/medium/split_$(s-1)_history.json"
			else
				jsonpath = "SubjectAnalysis_$crossOrWithin/deep/$t/medium/subject_$(s-1)_split_$(split)_history.json"
			end
			
			~isfile(jsonpath) ? continue : 1
			try
			d = JSON3.read(read(jsonpath))|>DataFrame
			d.split .= split
			d.subject .= s
			d.task .= t
			d.type .= crossOrWithin
			push!(dList,d)
			catch
				println("failed: ",t,s,split)
			end
		end
		
	end

# ╔═╡ c3a6bd53-da3b-470a-b723-b938fd1d6776
begin
tasklist = ["N170", "N400", "P3", "N2pc", "MMN", "ERN", "LRP"]
#tasklist = ["N170","P3"]#,"N400"]
dList = []
for t =tasklist
	for s =1:40
		
		read_json!(dList,t,s,"Cross")
		read_json!(dList,t,s,"Within")
	end
end
d = vcat(dList...)
	
d = select(d,Not(:batches))
end;

# ╔═╡ 732c10bc-4fd8-4f2c-9c61-d2c2b237df0a
md"""
we need to aggregate data across splits (for within subject)
"""

# ╔═╡ 1d0f4d1a-e144-45b9-9d25-c1c8c3d8aeff
d_agg = groupby(d,[:task,:subject,:epoch,:type]) |> 
		x->combine(x,
			:train_loss=>mean,
			:valid_loss=>mean,
			:valid_balanced_accuracy=>mean,
			renamecols=false);

# ╔═╡ c1e2557c-637a-4df9-9a23-b34cf30f6672
md"""
## Loading Encoding from erpcore
"""

# ╔═╡ fd48a6ff-8c35-426a-9f92-245558ecf4ba
begin
file = matopen("../data/D_amp_parent.mat")
erpcore = read(file,"D_p")
erpcore = .-diff(erpcore[:,:,:,1],dims=3)[:,:,1]
	
	erpcore = erpcore[:,:,1,1]
	
d_amp = DataFrame(erpcore,["P3","N170","MMN","N400","ERN","N2pc","LRP"])
d_amp.subject = 1:40
d_amp = stack(d_amp,variable_name="task",value_name="score")
sort!(d_amp,[:task,:subject])
end;

# ╔═╡ a4984b6d-6578-488b-aaa4-f4e46cd37e9f
data(d_amp) * mapping(:score,color=:task)*AlgebraOfGraphics.density()|>draw

# ╔═╡ 7d6bcbd6-e34d-4b5e-b50b-dba86e374aa6
md"""
#### Sidetrack amplitude vs difference
I am always looking for hints whether ERP amplitudes are additive or multiplicative. But I'd probably need to look at p2p vs. differences
"""

# ╔═╡ 7e8359fd-8c46-456c-9764-176959e98a84
let
	# plot mean amplitude vs difference, to check for potential multiplicative effects

	# => turns out, you should look at p2p for this to make sense.
	tmp = read(file,"D_p")


	d_diff = .-diff(tmp[:,:,:,1],dims=3)[:,:,1]
	d_mean = mean(tmp[:,:,:,1],dims=3)[:,:,1]
d_diff = DataFrame(d_diff,["P3","N170","MMN","N400","ERN","N2pc","LRP"])
d_diff.subject = 1:40
	d_diff.type .= "diff"
d_diff = stack(d_diff,variable_name="task",value_name="score")
	d_mean = DataFrame(d_mean,["P3","N170","MMN","N400","ERN","N2pc","LRP"])
d_mean.subject = 1:40
	d_mean.type .= "mean"
d_mean = stack(d_mean,variable_name="task",value_name="score")
	d = leftjoin(d_mean,d_diff,on=[:subject,:task],makeunique=true)
	#unstack(d,[],:type)
	@show d
	data(d)*mapping(:score=>"mean",:score_1=>"diff",color=:task,layout=:task)*visual(Scatter)|>x->draw(x,facet=(;linkxaxes=:none,linkyaxes=:none))
#scatter(tmp[:,:,1,1][:],tmp[:,:,2,1][:])


end

# ╔═╡ fd0d6584-63a8-45eb-8596-524c24823442
md"""
## load SME
these are uncertainty measures of the encoding amplitude/effects. I have this here because I accidentially loaded them thinking these were the encoding differences. I keep them around because they actually show some kind of correlation between tasks (e.g. noise is correlated across task), but due to space reason, I didnt include it in the paper
"""

# ╔═╡ c73d893d-5e41-4e35-9ff6-cf40f33ea624
begin
d_enc = []
for t = tasklist
	d = DataFrame(XLSX.readtable("../data/$(t)_Individual_SME.xlsx","Sheet1","J:K",first_row=4,header=true))
	d = d[1:end-4,:]
	d.task .= t
	push!(d_enc,d)
		
end
end

# ╔═╡ 308ace85-352e-4929-80cb-69137ab01bcd
md"""
## Combine loaded data to one dataframe
"""

# ╔═╡ f380980e-30bc-498a-a6b4-9926a196f813
begin
d_enc_agg = vcat(d_enc...)
d_dec = subset(d_agg,:epoch=>x->x.==maximum(d_agg.epoch));
#leftjoin(d_dec,d_enc;on=[:subject=>Symbol("Subject ID"),:task],matchmissing=:notequal)
rename!(d_enc_agg,Dict(Symbol("Subject ID")=>:subject,Symbol("Mean amplitude")=>:score))
sort!(d_enc_agg,[:task,:subject])
	
	d_enc_agg.type .= "sme"
	rename!(d_dec,Dict(:valid_balanced_accuracy => :score))
	#subset!(d_enc_agg,Not(:subject == 7 && :task=="N2pc")))
	#subset!(d_enc_agg,:score=>x->x.!="NaN")
	allowmissing!(d_enc_agg,:score)
	d_enc_agg.score[d_enc_agg.score .=="NaN"] .= missing
	d_fin = vcat(d_dec[:,[:subject,:task,:type,:score]],d_enc_agg)
	d_amp.type .= "enc"

	d_snr = deepcopy(d_amp)
	d_snr.score = d_amp.score ./ d_enc_agg.score
	d_snr.type .= "enc_snr"
	d_fin = vcat(d_fin,d_amp,d_snr)
	
	d_fin.score = convert.(Union{Float64,Missing},d_fin.score);
end;

# ╔═╡ 39f925d3-a8b8-42be-adc7-ce9f60387119
tlist = ["N170","P3","N400","MMN","N2pc","ERN","LRP"];
		

# ╔═╡ ab12ddc4-ff19-4d59-8789-763df8738e74
md"""
## Encoding vs. Decoding
"""

# ╔═╡ d48b3513-8c9b-4860-a876-1157472c09d4
begin
	f = Figure(resolution=(1000,400))
data(dropmissing(unstack(d_fin,:type,:score)))*mapping(:Cross=>(x->Float64.(x)) =>"balanced decoding accuracy [%]",:enc=>"Effect Amplitude [µV]",color=:task=>nonnumeric,layout=:task)*(visual(Scatter,alpha=0.5)+linear(;interval=nothing))|>x->draw!(f,x;axis=(;xticks=0.5:0.25:1),facet=(;linkyaxes=:none),palettes=(;layout=[(1,1),(1,2),(1,3),(2,1),(2,2),(2,3),(2,4)]))
	#ax = Axis(f[1,4])
	data(dropmissing(unstack(d_fin,:type,:score)))*mapping(:Cross=>(x->Float64.(x)) =>"",:enc=>"",color=:task=>nonnumeric)*(linear(;interval=nothing))|>x->draw!(f[1,4],x;axis=(;xticks=0.5:0.25:1,title="Summary"))
save("../plots/2023-02-24_DecVsEnc.svg",f)


	f
end

# ╔═╡ 2f4ef592-5e7a-4534-b074-e6dfb0cfd4de
md"""
## Calculate pairwise correlations
"""

# ╔═╡ 588ae9f5-907e-4937-b8a0-9de52692084b
md"""
## CorrelationMatrices
"""

# ╔═╡ dcf6096c-1c20-48c0-8190-f147bc563850
begin
d_unstack = unstack(d_fin,[:type,:subject],:task,:score)|>dropmissing
		transform!(d_unstack,AsTable(Not([:type,:subject]))=>ByRow(mean)=>:avg)
end;

# ╔═╡ c21d469a-e9e1-4f72-b4fe-b6eec4e67edd
begin
d_ρ = []

for (ix_col,t_col) = enumerate(tlist)
			for (ix_row,t_row) = enumerate(tlist)
				for typ = unique(d_unstack.type)
					if ix_col>=ix_row
						continue
					end
				d_loop = d_unstack[d_unstack.type.==typ,:]
				x = d_loop[:,t_col]
				y = d_loop[:,t_row]
				@rput x
				@rput y
				ρ = R"""
				library(WRS2)
				ρ = pbcor(x,y=y,ci=TRUE)
				"""
				#@show ρ["cor"]
				@rget ρ
				d_tmp = DataFrame(:corr=>ρ[:cor],:p_value=>ρ[:p_value],:cor_low=>ρ[:cor_ci][1],:cor_high=>ρ[:cor_ci][2],:task1=>t_col,:task2=>t_row,:type=>typ)
				push!(d_ρ,d_tmp)
				end
			end
end
d_ρ = vcat(d_ρ...)
end;


# ╔═╡ 1326254f-bdba-43c5-9153-b947eb0c37fe
data(d_ρ)*mapping(:corr,color=:type)*AlgebraOfGraphics.density()|>draw

# ╔═╡ c6da9a5f-ac5f-4761-8571-f0487c77c796
	let
		
		
		f = Figure()
		
		for (ix_col,t_col) = enumerate(tlist)
			for (ix_row,t_row) = enumerate(tlist)
				if ix_col==ix_row # switch to >= to drop diagonal :)
					# draw density
					d_tmp = unstack(subset(d_fin,:type=>x->x.!="Within",:task=>x->x.==t_col),:type,:score)
						grid = data(d_tmp)*mapping(:Cross,:enc,color=:subject)*visual(Scatter,markersize=7)|>x->draw!(f[ix_row,ix_row],x;axis=(ylabel="Decoding Accuracy",xticks=0.4:0.1:1,
					xticksize=4,
					xtickalign=1,))
					hidexdecorations!(current_axis(),ticks=false)
					hideydecorations!(current_axis(),ticks=true)
					
					text!(current_axis(),t_row,
						position=(0.05,0.8),align=(:left,:bottom),
						space=:relative,fontsize=10)
					
					continue
				elseif ix_col>ix_row
					typ = "enc"
					lim = nothing#Makie.Automatic
					
				else
					typ = "Cross"
					lim = (0.45,1.,0.45,1.)
				end
				d_loop = d_unstack[d_unstack.type.==typ,:]
				
				ρ = @subset( d_ρ,:type.==typ,:task1.==t_col,:task2.==t_row)
				if isempty(ρ) 
					ρ = @subset( d_ρ,:type.==typ,:task2.==t_col,:task1.==t_row)
				end
				pl = data(d_loop)*visual(Scatter,markersize=7)

				
				pl*mapping(t_col,t_row,color=:subject)|>
				x->AlgebraOfGraphics.draw!(
					f[ix_row,ix_col],x;
					axis=(
						#aspect=1,
						xticks=0.4:0.1:1,
						yticks=0.4:0.1:1,
					xticksize=4,
					yticksize=4,
					xtickalign=1,
					ytickalign=1,

						))
						
				isnothing(lim) ? "" : limits!(current_axis(),lim...) 
				text!(current_axis(),"$(round(ρ[1,:corr],digits=2)) $(round.(collect(ρ[1,[:cor_low,:cor_high]]),digits=2))",position=(0.1,0.9),space=:relative,fontsize=7)
				
				#text!(current_axis(),"bla",position=(0.7,0.7),space=:data)
				hidespines!(current_axis(), :t, :r) # only top and right
				if ix_col==1   && ix_row == length(tlist)
					hidexdecorations!(current_axis(),label=false,ticks=false)
					hideydecorations!(current_axis(),label=false,ticks=false)
				
				elseif ix_col==1  
					hidexdecorations!(current_axis())
					hideydecorations!(current_axis(),label=false,ticks=false)
				elseif ix_row == length(tlist)
					hidexdecorations!(current_axis(),label=false,ticks=false)
					hideydecorations!(current_axis())
				else
					hidedecorations!(current_axis())
				end
				#hidedecorations!(axright, grid = false)
				isnothing(lim) ? "" : ablines!(current_axis(),[0],[1];color="grey")
				#Axis
			colsize!(f.layout, ix_col, Aspect(1, 1))
			end
		end
		rowgap!(f.layout, Relative(0.02))
		colgap!(f.layout, Relative(0.02))
		resize_to_layout!(f)
		#CairoMakie.trim!(f.layout)
		#save("corrMat.svg", f)

		
		f
	end

# ╔═╡ 4f38c5ce-d96a-4962-a945-f4fae6bfec5c
let
	d = deepcopy(d_ρ)
	d.task1 = CategoricalArray(d.task1,levels=tlist)
	d.task2 = CategoricalArray(d.task2,levels=reverse(tlist))
	function clip_colormap(cm)
		alpha = abs.(range(-1,1,length=length(cm)))
		
		return [ColorTypes.RGBA.(c,alpha[i]>0.3 ? 1 : 0.3) for (i,c) in enumerate(cm)]
	end
	
data(@subset(d,:type.!=="Within",:type.!=="sme",:type.!=="enc_snr"))*mapping(:task1=>"",:task2=>"",:corr=>"robust ρ",layout=:type=>renamer("Cross"=>"Decoding","enc"=>"Encoding"))*
	visual(Heatmap,
		colorrange=(-1,1),
		colormap=reverse(ColorSchemes.get(ColorSchemes.RdBu,0:0.01:1))|>
				cm->clip_colormap(cm))|>
x->draw(x,axis=(;xgridvisible=false,ygridvisible=false,xticklabelrotation=π/2,aspect=1),figure=(;resolution=(500,250)))
	save("../plots/2023-02-04_corrHeatmap.svg", current_figure())
	current_figure()

end

# ╔═╡ 25b2002b-dbcf-43e5-9fa2-31f7ec3e7a74
IntervalsBetween

# ╔═╡ 9da73195-9ccb-493f-9a88-911dd4a977e3
begin
	mScore = @by d_fin [:type,:task] :meanScore=mean(abs.(skipmissing(:score)))
	
end;

# ╔═╡ fccf4782-83af-421c-af2b-e1f8fceed807

let
d_tmp = leftjoin(d_fin,mScore,on=[:type,:task])|>dropmissing
	d_tmp = d_tmp[d_tmp.type .=="Cross",:]
	d_tmp.aboveAvg = (d_tmp.score ./d_tmp.meanScore).>=1
	d_tmp2 =  @by d_tmp :subject :meanAboveAvg = mean(:aboveAvg)
	d_tmp = leftjoin(d_tmp,d_tmp2,on=:subject)
	sort!(d_tmp,:meanAboveAvg)
	
	d_tmp.subject = CategoricalArray(d_tmp.subject,levels=unique(d_tmp.subject))
	f = Figure(;resolution=(1000,250))
	ax =f[1,1] =  Axis(f;title="",xticksvisible=false,xticklabelsvisible=false,ylabel=L"\frac{\text{decoding accuracy}}{\text{mean taskwise accuracy}}",xlabel="subject, sorted according to overall decoding accuracy")
	hlines!(ax,[1],color="gray")
gh = data(d_tmp)*mapping(:subject=>levelcode,(:score,:meanScore)=>((x,y)->Float64.(x)./y))*(linear(interval=nothing)+mapping(color=:aboveAvg)*visual(Scatter))|>x->draw!(ax,x)
	AlgebraOfGraphics.legend!(f[1,2],gh)
	save("../plots/2023-02-27_aboveAVG.svg",current_figure())
	current_figure()
end

# ╔═╡ b6e1014a-5294-4043-9361-ed86581dc591


# ╔═╡ 7d399df0-6d61-4a17-9ca0-da78172a23d8
@by(d_ρ,:type,:max=maximum(:corr),:min=minimum(:corr),:maxp=maximum(:p_value),:minp=minimum(:p_value),:mean=mean(:corr))

# ╔═╡ 198a3056-aa64-443e-ba20-56308c7350b5
maximum(@subset(d_fin,:type.=="Cross").score)

# ╔═╡ 2393c535-c6e1-4d0c-b963-7ffb7153ec3b


# ╔═╡ 1b726dd3-e1b8-4cd9-807c-2ce2c783a20f
let
	d = deepcopy(d_ρ)
	d.task1 = CategoricalArray(d.task1,levels=tlist)
	d.task2 = CategoricalArray(d.task2,levels=reverse(tlist))
	
data(d)*mapping(:task1,:task2,:p_value=>"pvalue",col=:type)*visual(Heatmap,colorrange=(0.,0.1))|>x->draw(x,axis=(;xticklabelrotation=π/2,aspect=1))


end

# ╔═╡ edb8cf86-e1e3-4e75-8e73-a06538f589cd
begin
	gp = groupby(d_fin,:subject)
	ρ_singlesub = combine(gp;threads=false) do g
		dec = @subset(g,:type.=="Cross")
		enc = @subset(g,:type.=="enc")
		#ix = .!ismissing.(enc)
		#dec = dec[ix]
		#enc = enc[ix]
		decS = dec[in.(dec.task,Ref(enc.task)),:score]
		encS = enc[in.(enc.task,Ref(dec.task)),:score]
		
		@rput decS
		@rput encS
		ρ = R"""
				library(WRS2)
				ρ = pbcor(decS,y=encS,ci=TRUE)
				"""
		@rget ρ

		return ρ[:cor]
	end
#scatter(dec,enc)
end

# ╔═╡ 6dc5c553-7cfe-4b31-8189-5666bf5040bc
begin
	
	ρ_task = combine(groupby(d_fin,:task);threads=false) do g
		dec = @subset(g,:type.=="Cross")
		enc = @subset(g,:type.=="enc")
		#ix = .!ismissing.(enc)
		#dec = dec[ix]
		#enc = enc[ix]
		decS = dec[in.(dec.subject,Ref(enc.subject)),:score]
		encS = enc[in.(enc.subject,Ref(dec.subject)),:score]
		
		@rput decS
		@rput encS
		ρ = R"""
				library(WRS2)
				ρ = pbcor(decS,y=encS,ci=TRUE)
				"""
		@rget ρ
return DataFrame(:corr=>ρ[:cor],:p_value=>ρ[:p_value],:cor_low=>ρ[:cor_ci][1],:cor_high=>ρ[:cor_ci][2])
		
	end
#scatter(dec,enc)
end

# ╔═╡ 365a6a8f-c865-4d5e-8723-4409200a29fc
begin
hist(abs.(ρ_task.corr))
	vlines!([abs.(median(ρ_task.corr))])
	current_figure()
end

# ╔═╡ 1b32eb47-6782-4651-8c43-858e9fab8af6
	maximum(abs.(ρ_task.corr))

# ╔═╡ 2f28e2a8-d1d6-4b9c-abfe-2fc104380ed2


# ╔═╡ b2c01e22-91fb-4605-8792-3c7db68e26ff
begin
hist(ρ_singlesub.x1)
	vlines!([median(ρ_singlesub.x1)])
	current_figure()
end

# ╔═╡ 8e3202eb-de2f-42d9-9256-c133c032bc02


# ╔═╡ Cell order:
# ╠═6086b32e-a3ce-11ed-336c-1dd69a0fed42
# ╠═af862c8f-9d52-4d29-85fa-4289dae0c6d5
# ╟─4fde6495-cb2a-4855-af11-36ac865cd0d9
# ╠═19eb5619-c2d1-47d1-9c09-2484c9e9f00b
# ╠═c3a6bd53-da3b-470a-b723-b938fd1d6776
# ╟─732c10bc-4fd8-4f2c-9c61-d2c2b237df0a
# ╠═1d0f4d1a-e144-45b9-9d25-c1c8c3d8aeff
# ╠═c1e2557c-637a-4df9-9a23-b34cf30f6672
# ╠═fd48a6ff-8c35-426a-9f92-245558ecf4ba
# ╠═a4984b6d-6578-488b-aaa4-f4e46cd37e9f
# ╟─7d6bcbd6-e34d-4b5e-b50b-dba86e374aa6
# ╠═7e8359fd-8c46-456c-9764-176959e98a84
# ╟─fd0d6584-63a8-45eb-8596-524c24823442
# ╠═c73d893d-5e41-4e35-9ff6-cf40f33ea624
# ╟─308ace85-352e-4929-80cb-69137ab01bcd
# ╠═f380980e-30bc-498a-a6b4-9926a196f813
# ╠═39f925d3-a8b8-42be-adc7-ce9f60387119
# ╟─ab12ddc4-ff19-4d59-8789-763df8738e74
# ╠═d48b3513-8c9b-4860-a876-1157472c09d4
# ╟─2f4ef592-5e7a-4534-b074-e6dfb0cfd4de
# ╠═c21d469a-e9e1-4f72-b4fe-b6eec4e67edd
# ╠═1326254f-bdba-43c5-9153-b947eb0c37fe
# ╠═c6da9a5f-ac5f-4761-8571-f0487c77c796
# ╟─588ae9f5-907e-4937-b8a0-9de52692084b
# ╠═dcf6096c-1c20-48c0-8190-f147bc563850
# ╠═4f38c5ce-d96a-4962-a945-f4fae6bfec5c
# ╠═25b2002b-dbcf-43e5-9fa2-31f7ec3e7a74
# ╠═9da73195-9ccb-493f-9a88-911dd4a977e3
# ╠═fccf4782-83af-421c-af2b-e1f8fceed807
# ╠═b6e1014a-5294-4043-9361-ed86581dc591
# ╠═7d399df0-6d61-4a17-9ca0-da78172a23d8
# ╠═198a3056-aa64-443e-ba20-56308c7350b5
# ╠═2393c535-c6e1-4d0c-b963-7ffb7153ec3b
# ╠═1b726dd3-e1b8-4cd9-807c-2ce2c783a20f
# ╠═edb8cf86-e1e3-4e75-8e73-a06538f589cd
# ╠═6dc5c553-7cfe-4b31-8189-5666bf5040bc
# ╠═365a6a8f-c865-4d5e-8723-4409200a29fc
# ╠═1b32eb47-6782-4651-8c43-858e9fab8af6
# ╠═2f28e2a8-d1d6-4b9c-abfe-2fc104380ed2
# ╠═b2c01e22-91fb-4605-8792-3c7db68e26ff
# ╠═8e3202eb-de2f-42d9-9256-c133c032bc02
