import nmrglue as ng
import matplotlib
import numpy as np
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

dic, data = ng.bruker.read('/Users/mnk/Downloads/NMR-2/jn290-1A/1')
data = ng.bruker.remove_digital_filter(dic, data)


# process the spectrum
data = ng.proc_base.zf_size(data, 32768)    # zero fill to 32768 points
data = ng.proc_base.fft(data)               # Fourier transform
data = ng.proc_autophase.autops(data, 'peak_minima')     # phase correction
data = ng.proc_base.di(data)                # discard the imaginaries
data = ng.proc_base.rev(data)
# reverse the data



udic = ng.bruker.guess_udic(dic, data)
uc = ng.fileiobase.uc_from_udic(udic)
ppm_scale = uc.ppm_scale()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(ppm_scale, data)

# read in the integration limits
peak_list = open("/Users/mnk/Downloads/NMR-2/jn290-1A/1/pdata/1/intrng", "r")

peak_list.readline()
for line in peak_list:
    x = line.split("   ")

    start = float(x[0])
    end = float(x[1][:-1])
    min = uc(start, "ppm")
    max = uc(end, "ppm")
    if min > max:
        min, max = max, min

    # extract the peak
    peak = data[min:max + 1]
    peak_scale = ppm_scale[min:max + 1]

    # plot the integration lines, limits and name of peaks
    ax.plot(peak_scale, peak.cumsum() / 100. + peak.max())
    ax.plot(peak_scale, [0] * len(peak_scale))
    ax.text(peak_scale[0], 0.5 * peak.sum() / 100. + peak.max(), peak.sum(),fontsize=8)


peak_list.close()
ax.set_xlim(ppm_scale[0], ppm_scale[-1])
fig.savefig('figure_nmrglue.png')