import nmrglue as ng
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import functions as f

path = "/Users/mnk/Downloads/NMR-2/jn290-1B/1/pdata/1"
impurities = {"CDCl3": {"solvent": 7.26, "H2O": 1.56}}

dic, data = ng.bruker.read("/Users/mnk/Downloads/NMR-2/jn290-1B/1")
data = ng.bruker.remove_digital_filter(dic, data)

# process the spectrum
data = ng.proc_base.zf_size(data, 32768)  # zero fill to 32768 points
data = ng.proc_base.fft(data)  # Fourier transform
data = ng.proc_autophase.autops(data, 'peak_minima')  # phase correction
data = ng.proc_base.di(data)  # discard the imaginaries
data = ng.proc_base.rev(data)
# reverse the data

peak_table = ng.peakpick.pick(data, pthres=1e6, algorithm='downward')

# conversion to ppm
udic = ng.bruker.guess_udic(dic, data)
uc = ng.fileiobase.uc_from_udic(udic)
ppm_scale = uc.ppm_scale()


fig = plt.figure()
ax = fig.add_subplot(111)

# read in the integration limits
peak_list = open(path + "/intrng", "r")

solvent_file = open(path + "/proc", "r");
solvent = ""
for line in solvent_file:
    if "##$SREGLST= " in line:
        solvent = line[16:-2]

integral_list = {}
peak_list.readline()  # delete first unecessarry line
cont = False

for line in peak_list:
    x = line.split("   ")
    b = float(x[1][:-1])
    e = float(x[0])

    start = uc(b, "ppm")
    end = uc(e, "ppm")
    if start > end:
        start, end = end, start

    for i in impurities[solvent]:
        if (b < impurities[solvent][i] < e) or (b < 0 or e < 0):
            data[start:end] = 0
            cont = True

    if cont:
        cont = False
        continue

    # extract the peak
    peak = data[start:end + 1]
    integral_list[float(peak.sum())] = [float(start), float(end), []]
peak_list.close()


peaks = ng.peakpick.pick(data, 1e10)
print(len(peaks))
for peak in peaks['X_AXIS']:
    for i in integral_list:
        if integral_list[i][0] < peak < integral_list[i][1]:
            integral_list[i][2].append(peak)

#filter out 5%
percent = 5
integral_list = f.delete_impurities(integral_list, percent, data)

#ratios of integrals
integral_list = f.find_ratios(integral_list)

#draw integrals
f.draw_integrals(integral_list, data, ppm_scale, ax)

#draww plot
ax.plot(ppm_scale, data, color="black")
ax.set_xlim(ppm_scale[0], ppm_scale[-1])
fig.savefig('figure_nmrglue.png')


print(integral_list)
peak_locations_ppm = [uc.ppm(i) for i in peaks['X_AXIS']]
peak_amplitudes = data[peak_table['X_AXIS'].astype('int')]


