# PANDORAS

CMAQ/PANDORAS evaluation system

## Description and Notes

Extracts CMAQ at PANDORAS locations and times. PANDORAS locations include
both the base location and the offset due to light path slant. CMAQ output
is the column density in dobson units.

Assumptions:

1. Only using SZA<70 degrees,
2. Simple cartesian coordinates are aligned with the model grid,
3. Vertical structure is constant within slant-based horizontal offsets,
4. Horizontal cell at ZH is representative from ZF[l-1] to ZF[l],

Limitations:

1. I currently am not doing any treatment of the stratosphere,
2. This isn't production, so I am not doing any error checking.
3. Currently, only the model, observation and observation uncertainty are output


Findings using a 4km grid:

1. Column could be displaced by 12 and Row by 5 within CMAQ  vertical grid
2. Column could be displaced by 3 and Row by 1 within first 4km
3. Results depend on SZA limits (see Assumptions)

## Prerequisites

* PANDORAS L2 file
* CMAQ conc file
* METCRO3D file (or equivalent with PRES and TA)
* Python (>=3.6) with PseudoNetCDF

## Commands To Reproduce

Edit Makefile as appropriate
```
make
```

## Directory Structure

figs and output contents only exist after successful make

```
.
|-- README.md
|-- Makefile
|-- output
|   `-- SchillerParkIL_%Y%m%d.nc
|-- figs
|   |-- SchillerParkIL_20170522-20170618.hourly.png
|   |-- SchillerParkIL_20170522-20170618.png
|   |-- SchillerParkIL_%Y%m%d.hourly.png
|   `-- SchillerParkIL_%Y%m%d.png
|-- lb3.pandonia.net
|   `-- SchillerParkIL
|       `-- Pandora51s1
|           `-- L2
|               `-- Pandora51s1_SchillerParkIL_L2Tot_rnvs0p1-5.txt
`-- scripts
    |-- pair.py
    |-- pandoras.py
    `-- plot.py
```
