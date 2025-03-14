
# Debug script for peakfiles. 
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendDSP
using LegendDSP: get_fltpars
using LegendSpecFits: get_friedman_diaconis_bin_width
using LegendHDF5IO
using RadiationSpectra
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using LegendMakie, Makie, CairoMakie
using Measures

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/processing_funcs/process_peak_split.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/src/apply_qc.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

# inputs
reprocess = true
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(2)
channel = ChannelId(1)
category = DataCategory(:cal)
source = :co60

# load configs 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
qc_config = dataprod_config(asic).qc(filekeys[1]).default
ecal_config = dataprod_config(asic).energy(filekeys[1]).default
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)

# DEBUG
data = asic
filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])

if Symbol(ecal_config.source) == :co60
    gamma_lines =  sort(ecal_config.co60_lines)
    gamma_names =  ecal_config.co60_names
    left_window_sizes = ecal_config.co60_left_window_sizes 
    right_window_sizes = ecal_config.co60_right_window_sizes 
elseif Symbol(ecal_config.source) == :th228
    gamma_lines =  ecal_config.th228_lines
    gamma_names =  ecal_config.th228_names
    left_window_sizes = ecal_config.left_window_sizes 
    right_window_sizes = ecal_config.right_window_sizes 
end

h_uncals = Vector{Histogram}(undef, length(filekeys))
peakpos = Vector{Vector{<:Real}}(undef, length(filekeys))
f = 1
# for f in eachindex(filekeys) 
filekey = filekeys[f]
# peak_file = peak_files[f]
data_ch = read_ldata(data, DataTier(:raw), filekey, channel)
wvfs = data_ch.waveform
e_uncal = filter(x -> x >= qc_config.e_trap.min , data_ch.daqenergy)
if isempty(e_uncal)
    @warn "No energy values >= $(qc_config.e_trap.min) found for $filekey - skip"
    # continue
end
eventnumber = data_ch.eventnumber
timestamp = data_ch.timestamp

    # binning and peak search windows and histogram settings 
bin_min = quantile(e_uncal, ecal_config.left_bin_quantile)
bin_max = quantile(e_uncal, ecal_config.right_bin_quantile)
peak_min = quantile(e_uncal, ecal_config.left_peak_quantile)
peak_max = quantile(e_uncal, ecal_config.right_peak_quantile)
bin_width = LegendSpecFits.get_friedman_diaconis_bin_width(filter(in(bin_min..bin_max), e_uncal))
if (peak_max-peak_min)/bin_width < ecal_config.nbins_min
    bin_width = (peak_max-peak_min)/ ecal_config.nbins_min
end

# peak search
h_uncals[f] = fit(Histogram, e_uncal, 0:bin_width:maximum(e_uncal)) # histogram over full energy range; stored for plot 

try
    h_peaksearch = fit(Histogram, e_uncal, bin_min:bin_width:peak_max) # histogram for peak search
    _, peakpos[f] = RadiationSpectra.peakfinder(h_peaksearch, σ= ecal_config.peakfinder_σ, backgroundRemove=true, threshold = ecal_config.peakfinder_threshold)
catch e 
    @warn "peakfinder failed for $filekey - use larger window. err $e"
    h_peaksearch = fit(Histogram, e_uncal, bin_min:bin_width:(peak_max*1.5)) # histogram for peak search
    _, peakpos[f] = RadiationSpectra.peakfinder(h_peaksearch, σ= ecal_config.peakfinder_σ, backgroundRemove=true, threshold = ecal_config.peakfinder_threshold)
end 
if length(peakpos[f]) !== length(gamma_lines)
    error("Number of peaks found $(length(peakpos[f])); expected gamma lines $(length(gamma_lines)) \n you could try to modify peakfinder_threshold and/or peakfinder_σ (file $f)")
else 
    @info "Found $(length(peakpos[f])) peaks for $filekey"
end 
cal_simple = mean(gamma_lines./sort(peakpos[f]))
e_simplecal = e_uncal .* cal_simple
i = 1
peakIdx = findall((gamma_lines[i] - left_window_sizes[i]) .<= e_simplecal .< (gamma_lines[i] + right_window_sizes[i]))
# do simple dsp. only for peaks. 
dsp_par = simple_dsp_qc(Table(waveform = wvfs[peakIdx]), dsp_config)

