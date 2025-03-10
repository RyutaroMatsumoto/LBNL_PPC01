using LegendDataManagement
using LegendDataManagement.LDMUtils
using CairoMakie, Plots
using LegendHDF5IO
using Unitful
using TypedTables
using StatsBase
using ColorSchemes

asic = LegendData(:ppc01)
period = DataPeriod(3)
category = :cal 
channel = ChannelId(1)


runs1 = [18, 19, 20, 21]
runs2 = [23, 24, 25, 26]

taw = [0.04, 0.08, 0.16, 0.32]
fwhm_Co60a1 = [asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fit.Co60a.fwhm for r in runs1]
fwhm_Co60a2 = [asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fit.Co60a.fwhm for r in runs2]
fwhm_Co60b1 = [asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fit.Co60b.fwhm for r in runs1]
fwhm_Co60b2 = [asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fit.Co60b.fwhm for r in runs2]
peaks =  round.(ustrip.(asic.par.rpars.ecal[period, DataRun(21), channel].e_trap.cal.peaks), digits = 1 )
palette = ColorSchemes.twelvebitrainbow
colors = [get(palette, i/length(runs)) for i in eachindex(runs)]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Trigger Averaging Window (µs)", ylabel = "FWHM (keV)", xscale = log10, xticks = (taw, string.(taw)), #(3.4, 5.3)
                title = "LBNL-PPC01 and ASIC L1k65n with buffer:  $period, r00$(runs1[1])-r0$(runs1[end])", titlesize = 16)
for i in eachindex(runs1)
    CairoMakie.scatter!(ax, taw[i], ustrip.(fwhm_Co60a1[i]), color = colors[1], label = "$(Makie.UnicodeFun.to_superscript(60))Co a) $(peaks[1]) keV")
    CairoMakie.scatter!(ax, taw[i], ustrip.(fwhm_Co60b1[i]), color = colors[end-1], label = "$(Makie.UnicodeFun.to_superscript(60))Co b) $(peaks[2]) keV")
end 
axislegend("Gamma peak", position = :lb, merge = true)
fig
plt_folder = "/global/homes/r/ryutaro/lbnl/PPC01_ryutaro/generated/jlplt/rplt/cal/p03/compare"
if !isdir(plt_folder)
    mkdir(plt_folder)
end
plt_name = joinpath(plt_folder, "taw_fwhm_8.png")
save(plt_name, fig)

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Pulse Height Averaging Window (µs)", ylabel = "FWHM (keV)", xscale = log10, xticks = (taw, string.(taw)), #(3.4, 5.3)
                title = "LBNL-PPC01 and ASIC L1k65n with buffer:  $period, r00$(runs2[1])-r0$(runs2[end])", titlesize = 16)
for i in eachindex(runs2)
    CairoMakie.scatter!(ax, taw[i], ustrip.(fwhm_Co60a2[i]), color = colors[1], label = "$(Makie.UnicodeFun.to_superscript(60))Co a) $(peaks[1]) keV")
    CairoMakie.scatter!(ax, taw[i], ustrip.(fwhm_Co60b2[i]), color = colors[end-1], label = "$(Makie.UnicodeFun.to_superscript(60))Co b) $(peaks[2]) keV")
end 
axislegend("Gamma peak", position = :rt, merge = true)
fig

plt_name = joinpath(plt_folder, "taw_fwhm_32.png")
save(plt_name, fig)