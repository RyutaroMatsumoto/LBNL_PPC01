using LegendDataManagement
using Makie, CairoMakie, LegendMakie
using Measurements: value as mvalue, uncertainty as muncert
using Unitful 
using StatsBase

# setup 
asic = LegendData(:ppc01)
period = DataPeriod(3)
channel = ChannelId(1)
category = :cal 
e_types = [:e_trap]
runs = 50:52

# load pars 
cable_length = sort([185.0, 405.0, 645.0 ]; rev = false)
peak = :Tl208FEP
fwhm_peak = ustrip.([asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fit[peak].fwhm for r in runs])
fwhm_qbb = ustrip.([asic.par.rpars.ecal[period, DataRun(r), channel].e_trap.fwhm.qbb for r in runs])

# plot 
fig = Figure()
ax = Axis(fig[1,1], xticks = cable_length,  
                    xlabel = "Cable length (cm)", ylabel = "FWHM (keV)")
hlines!(ax, [mvalue(mean(fwhm_peak))], linestyle = :dash, alpha = 0.5)
errorbars!(ax, cable_length, mvalue.(fwhm_peak), muncert.(fwhm_peak), label = "Tl208FEP (2.614 MeV)", whiskerwidth = 0) 
scatter!(ax, cable_length, mvalue.(fwhm_peak), label = "Tl208FEP (2.614 MeV)", markersize = 10) 
hlines!(ax, [mvalue(mean(fwhm_qbb))], linestyle = :dash, alpha = 0.5)
errorbars!(ax, cable_length, ustrip.(mvalue.(fwhm_qbb)), muncert.(fwhm_qbb), label = "Qbb (2.039 MeV)", whiskerwidth = 0) 
scatter!(ax, cable_length, mvalue.(fwhm_qbb), label = "Qbb (2.039 MeV)", markersize = 10) 
axislegend(merge = true, orientation = :horizontal, position = :lt)

Makie.ylims!(ax, 3.4 , 4.6 )
fig
save("$(@__DIR__)/plots/Th228_resolution_vs_cablelength.png", fig)

