# -*- coding: latin1 -*-
import PseudoNetCDF as pnc
from collections import OrderedDict
import pandas as pd
import numpy as np
import re


"""
_dates = {}


def pdate(x):
    if x not in _dates:
        _dates[x] = datetime.strptime(x + '+0000', '%Y%m%dT%H%M%SZ%z')

    return _dates[x]
"""

_keyval = re.compile(r'(.+?)\s*:\s*(.+)')
_units = re.compile(r'.*\[(.+)\].*')
_col2name = {
    'Column 1': 'time',
    'Column 2': 'FDAY',
    'Column 3': 'DURATION',
    'Column 4': 'SolarZenithAngle',
    'Column 5': 'SolarAzimuthAngle',
    'Column 6': 'LunarZenithAngle',
    'Column 7': 'LunarAzimuthAngle',
    'Column 8': 'SpeciesColumnAmount',
    'Column 9': 'SpeciesColumnUncertainty',
    'Column 10': 'AirMassFactor',
    'Column 11': 'DiffuseCorrection',
    'Column 12': 'QualityAssuranceFlag',
    'Column 13': 'DQ1Flags',
    'Column 14': 'DQ2Flags',
    'Column 15': 'FittingResultIndex',
    'Column 16': 'NormalizedRMS',
    'Column 17': 'ExpectedNormalizedRMSMeasurement',
    'Column 18': 'ExpectedNormalizedRMSInstrument',
    'Column 19': 'ClimatologicalStationPressure',
    'Column 20': 'DataProcessingTypeIndex',
    'Column 21': 'CalibrationFileVersion',
    'Column 22': 'CalibrationFileValidityStartingDate',
    'Column 23': 'Level2FitDataQualityFlag',
    'Column 24': 'L2DQ1Flag',
    'Column 25': 'L2DQ2Flag',
    'Column 26': 'Level1DataQualityFlag',
    'Column 27': 'L1DQ1Flag',
    'Column 28': 'L1DQ2Flag',
    'Column 29': 'WavelengthEffectiveTemperature',
    'Column 30': 'EstimatedAverageResidualStrayLightLevel',
    'Column 31': 'RetrievedWavelengthShiftL1',
    'Column 32': 'RetrievedWavelengthShiftSpectral',
    'Column 33': 'IntegrationTime',
    'Column 34': 'NumberDarkCountCycles',
    'Column 35': 'EffectivePositionFilterwheel1',
    'Column 36': 'EffectivePositionFilterwheel2',
}


class pandoraslb3(pnc.PseudoNetCDFFile):
    def __init__(self, path):
        """
        Arguments
        ---------
        path : str
            path to pandoras file of lb3 format (see lb3.pandonia.net)

        Returns
        -------
        None
        """
        infile = open(path, mode='r', encoding='latin1')
        props = OrderedDict()
        skiplines = 0
        while True:
            skiplines += 1
            line = infile.readline()
            if line.startswith('-------'):
                break
            key, val = _keyval.match(line).groups()
            try:
                val = eval(val)
            except Exception:
                pass
            props[key] = val
        varprops = OrderedDict()
        while True:
            skiplines += 1
            line = infile.readline()
            if line.startswith('-------'):
                break
            key, val = _keyval.match(line).groups()
            key = key.strip()
            val = val.strip()
            varprops[_col2name.get(key, key)] = val

        varkeys = list(varprops)
        data = pd.read_csv(
            path, skiprows=skiplines, names=varkeys,
            parse_dates=['time'], delimiter=' '
        )
        self.setncatts(props)
        self.createDimension('time', data.shape[0])
        for varkey, varprop in varprops.items():
            var = self.createVariable(varkey, 'd', ('time',))
            var.description = varprop
            m = _units.match(varprop)
            if m is not None:
                var.units, = m.groups()
            var[:] = data[varkey]

        self.variables['time'][:] /= 1e9
        self.variables['time'].units = (
            'seconds since 1970-01-01 00:00:00+0000'
        )

    def ray2txyz(self, toa=20000, zincr=None):
        """
        Description
        -----------
        Convert time, SolarZenithAngle and SolarAzimuthAngle to
        coordinates t (seconds since epoch)

        Arguments
        ---------
        toa : int
            top of the atmosphere in meters; ray segments will be long
            enough so that the returned z coordinate reaches toa
        zincr : int (optional)
            if zincr is provided, calculate from 0m to max legnth of ray
            segment in incrmements of zincr

        Returns:
            t, x, y, z : array of coordinates in seconds (t) and meters
        """
        theta = np.radians(
            np.ma.masked_less(
                np.ma.masked_greater(
                    self.variables['SolarZenithAngle'][:],
                    70
                ),
                -70
            )
        )
        phi = np.ma.masked_array(
            np.radians(self.variables['SolarAzimuthAngle'][:])
        )
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        cosphi = np.cos(phi)
        sinphi = np.sin(phi)
        rtop = toa / costheta
        if zincr is None:
            r = np.ma.array([0 * costheta, rtop])
        else:
            r = np.arange(0, rtop.max() + zincr, zincr)[:, None]

        # r = np.arange(0, zmax, zincr)
        # https://en.wikipedia.org/wiki/Horizontal_coordinate_system
        # https://en.wikipedia.org/wiki/Spherical_coordinate_system
        # Switching assignment of x and y  from typical spherical coordinate
        # because Azimuth is measured 0 at N and latitude is measured 90 at N
        y = r[:] * sintheta * cosphi
        x = r[:] * sintheta * sinphi
        z = r[:] * costheta
        t = self.variables['time'][None, :].repeat(r.shape[0], 0)

        return t, x, y, z

    def findcells(self, metcro3df):
        """
        Arguments
        ---------
        metcro3df : PseudoNetCDF-like
            file must contain ZH variable and provide time2t, ll2ij, and ll2xy

        Returns
        -------
        t, k, j, i : array of coordinates to sample IOAPI file

        Note Assumptions:
            1. ZF(t, l, j, i) ~ ZF(t, l, j\pm 3, i\pm 3)
            2. TERRAIN(j, i) ~ TERRAIN(j\pm 3, i\pm 3)
            3. ZF centers represent slope well enough
        """
        t, dxs, dys, dzs = self.ray2txyz()
        mylat = getattr(self, 'Location latitude [deg]')
        mylon = getattr(self, 'Location longitude [deg]')
        ts = self.getTimes()
        tis = metcro3df.time2t(ts)
        i0, j0 = metcro3df.ll2ij(mylon, mylat)
        x0, y0 = metcro3df.ll2xy(mylon, mylat)
        print(i0, j0)
        ZH = metcro3df.variables['ZH'][:, :, j0, i0]
        coords = []
        k = np.arange(ZH.shape[1], dtype='i')
        for ti, dx, dy, dz in zip(tis, dxs.T, dys.T, dzs.T):
            zh = ZH[ti]
            ldx = np.interp(zh[:], dz, dx)
            ldy = np.interp(zh[:], dz, dy)
            x = x0 + ldx
            y = y0 + ldy
            i, j = metcro3df.ll2ij(*metcro3df.xy2ll(x, y))
            coords.append((ti * np.ones_like(k), k, j, i))
        t, k, j, i = np.array(coords, dtype='i').swapaxes(0, 1)
        return t, k, j, i
