import pandoras
import PseudoNetCDF as pnc
from datetime import datetime
import sys
import os
import gc


exists = os.path.exists

# dummy assignent
f = pandoras
panpath = sys.argv[1]
concpat = sys.argv[2]
metpat = sys.argv[3]
outpat = sys.argv[4]


allf = pnc.pncopen(panpath, format='pandoraslb3')
allptimes = allf.getTimes()
moddates = sorted([
    datetime.strptime(p, metpat)
    for p in list(set([
        t.strftime(metpat) for t in allptimes
    ]))
])
for moddate in moddates:
    metpath = moddate.strftime(metpat)
    concpath = moddate.strftime(concpat)
    outpath = moddate.strftime(outpat)
    if not exists(metpath):
        print('Skipping missing met date', metpath, flush=True)
        continue
    elif not exists(concpath):
        print('Skipping missing model date', concpath, flush=True)
        continue
    elif exists(outpath):
        print('Keeping cached output', outpath, flush=True)
        continue
    else:
        print('Processing', moddate.strftime('%F'), flush=True)

    mf = pnc.pncopen(metpath, format='ioapi').subset(
        ['ZF', 'PRES', 'TA', 'ZH']
    ).eval("""
DZ = ZF * 1
DZ[:, 1:] = np.diff(ZF[:], axis=1)
PPM2DU = 1 / 1e6 * PRES[:] / 8.314 / TA[:] * DZ[:] / 1e4 * 6.022e23 / 2.69e16
ZH = ZH
ZF = ZF
""")
    cf = pnc.pncopen(concpath, format='ioapi').subset(
        ['NO2']
    )
    cf.copyVariable(mf.variables['PPM2DU'], key='PPM2DU')
    cf.copyVariable(mf.variables['ZH'], key='ZH')
    cf.eval(
        """
NO2 = NO2[:] * PPM2DU
NO2.units = 'Dobsons'
NO2.var_desc = 'NO2 columns'
""",
        inplace=True
    )
    modts = cf.getTimes()
    stime = modts[0]
    etime = modts[-1]
    f = allf.slice(time=(allptimes >= stime) & (allptimes < etime))
    t, k, j, i = f.findcells(cf)
    ts = f.getTimes()
    mod = cf.variables['NO2'][t, k, j, i].sum(1)
    mod.dimensions = ('time',)
    sza = f.variables['SolarZenithAngle']
    saa = f.variables['SolarAzimuthAngle']
    obs = f.variables['SpeciesColumnAmount']
    obsu = f.variables['SpeciesColumnUncertainty']
    outf = pnc.PseudoNetCDFFile()
    outf.createDimension('time', len(ts))
    outf.createDimension('layer', len(cf.dimensions['LAY']))
    outf.copyVariable(f.variables['time'], key='time')
    vari = outf.createVariable('i', 'f', ('time', 'layer'))
    vari.long_name = vari.var_desc = 'i'
    vari.units = '0-based index'
    vari[:] = i
    varj = outf.copyVariable(vari[:] * 0 + j, key='j')
    varj.long_name = varj.var_desc = 'j'
    vark = outf.copyVariable(vari[:] * 0 + k, key='k')
    vark.long_name = vark.var_desc = 'k'
    outf.copyVariable(saa, key='SAA')
    outf.copyVariable(sza, key='SZA')
    outf.copyVariable(mod, key='MOD')
    outf.copyVariable(obs, key='OBS')
    outf.copyVariable(obsu, key='OBSU')
    outf.save(outpath)
    gc.collect()
