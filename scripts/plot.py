# from matplotlib import use; use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import PseudoNetCDF as pnc
import pandas as pd
import sys


figpath = sys.argv[1]
inpaths = sys.argv[2:]

inf = pnc.pncmfopen(inpaths, format='netcdf', stackdim='time')
ts = inf.getTimes()
obs = inf.variables['OBS'][:]
obsu = inf.variables['OBSU'][:]
mod = inf.variables['MOD'][:]
lb = obs - 1.96 * np.abs(obsu)
tsm = np.ma.masked_where(lb < 0, ts).compressed()
obsm = np.ma.masked_where(lb < 0, obs).compressed()
obsum = np.ma.masked_where(lb < 0, obsu).compressed()
modm = np.ma.masked_where(lb < 0, mod).compressed()
df = pd.DataFrame(dict(obs=obsm, obsu=obsum, mod=modm), index=tsm).resample('H').mean()
lr = scipy.stats.mstats.linregress(
    np.ma.masked_invalid(df['mod'].values),
    np.ma.masked_invalid(df['obs'].values),
)
ax = df.plot(y=['obs', 'mod'], linestyle='none', marker='o')
ax.text(0, 1, 'r={:.2f}'.format(lr.rvalue), transform=ax.transAxes)
ax.figure.savefig(figpath.replace('.png', '.hourly.png'))
plt.close()

fig = plt.figure(figsize=(16, 4))
ax = fig.add_axes([.1, .2, .85, .75])

ax.set_title('Schiller Park NO2 Columns')
ax.errorbar(tsm, obsm, 1.96 * obsum, color='k', linestyle='none', marker='o', zorder=1)
#ax.plot(
#    tsm, obsm, color='k', linestyle='none', marker='o',
#    label='obs {:.2g}$\pm${:.2g}'.format(float(obs.mean()), float(obs.std()))
#)
lr = scipy.stats.mstats.linregress(modm, obsm)
ax.plot(
    tsm, modm, linestyle='none', marker='o', color='r', zorder=3,
    label='mod {:.2g}$\pm${:.2g} (r={:.2f})'.format(float(mod.mean()), float(mod.std()), lr.rvalue)
)
#ax.set_ylim(0.05, 10)
# ax.set_yscale('log')
plt.legend()
plt.setp(ax.get_xticklabels(), rotation=45)
plt.savefig(figpath)
