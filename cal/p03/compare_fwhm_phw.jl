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
runs = [34,35,36]

phw = [5.00, 6.00, 4.00]
fwhm_Co60a = [asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fit.Co60a.fwhm for r in runs]
fwhm_Co60b = [asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fit.Co60b.fwhm for r in runs]
peaks =  round.(ustrip.(asic.par.rpars.ecal[period, DataRun(3), channel].e_trap.cal.peaks), digits = 1 )
palette = ColorSchemes.twelvebitrainbow
colors = [get(palette, i/length(runs)) for i in eachindex(runs)]

fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Pulse Height Window (µs)", ylabel = "FWHM (keV)", xscale = log10, xticks = (phw, string.(phw)), #(3.4, 5.3)
                title = "LBNL-PPC01 and ASIC L1k65n with buffer:  $period, r00$(runs[1])-r0$(runs[end])", titlesize = 16)
for i in eachindex(runs)
    CairoMakie.scatter!(ax, phw[i], ustrip.(fwhm_Co60a[i]), color = colors[1], label = "$(Makie.UnicodeFun.to_superscript(60))Co a) $(peaks[1]) keV")
    CairoMakie.scatter!(ax, phw[i], ustrip.(fwhm_Co60b[i]), color = colors[end-1], label = "$(Makie.UnicodeFun.to_superscript(60))Co b) $(peaks[2]) keV")
end 
axislegend("Gamma peak", position = :rc, merge = true)
fig

plt_folder = "/global/homes/r/ryutaro/lbnl/PPC01_ryutaro/generated/jlplt/rplt/cal/p03/compare"
if !isdir(plt_folder)
    mkdir(plt_folder)
end
plt_name = joinpath(plt_folder, "phw_fwhm_2m_r034-r036.png")
save(plt_name, fig)


