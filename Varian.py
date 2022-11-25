import nmrglue as ng
import matplotlib
import numpy as np
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


dir_name = '/Users/mnk/Downloads/NMR-2/jn197-3B_20210212_01/PROTON_01.fid' # directory where the Varian data is held.
dic, data = ng.varian.read(dir_name, procpar_file='procpar')

udic = ng.varian.guess_udic(dic, data)
udic[0]['size']     = 1500
udic[0]['complex']  = True
udic[0]['encoding'] = 'direct'
udic[0]['sw']       = 50000.0
udic[0]['obs']      = 125.681
udic[0]['car']      = 99.0*125.681
udic[0]['label']    = '1H'

C = ng.convert.converter()
C.from_varian(dic, data, udic)
pdic, pdata = C.to_pipe()

pdic, pdata = ng.pipe_proc.sp(pdic, pdata, off=0.35, end=0.98, pow=2, c=1.0)
pdic, pdata = ng.pipe_proc.zf(pdic, pdata, auto=True)
pdic, pdata = ng.pipe_proc.ft(pdic, pdata, auto=True)
pdata = ng.proc_autophase.autops(pdata, 'peak_minima')
pdic, pdata = ng.pipe_proc.di(pdic, pdata)

# create a unit conversion object for the axis
uc = ng.pipe.make_uc(pdic, pdata)

# plot the spectrum
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(uc.ppm_scale(), pdata, 'k-')
#ax.set_xlim(200, -100)

fig.savefig('figure_nmrglue.png')