# L2FILE is the base PANDORAS file
L2FILE=lb3.pandonia.net/SchillerParkIL/Pandora51s1/L2/Pandora51s1_SchillerParkIL_L2Tot_rnvs0p1-5.txt

# CONC and METCRO3D are CMAQ concentration and met files
CONCTMPL=/work/REGIONS/users/jliljegr/lmos_2017/4LMOS3.baseline2017/CMAQv521_cb6r3_ae6nvPOA_aq.4LMOS3.35.4LMOS3.baseline2017.CONC.%Y%j
METCRO3DTMPL=/work/ROMO/finescale/cmaq/4LMOS3/met/METCRO3D.4LMOS3.35L.%y%m%d

# OUTROOT is the base file name for outputs and figures
OUTROOT=SchillerParkIL
OUTTMPL=output/$(OUTROOT)_%Y%m%d.nc
FIGIN=output/$(OUTROOT)_????????.nc
FIGOUT=figs/$(OUTROOT)_20170522-20170618.png

all: $(L2FILE) output/updated $(FIGOUT)

output/updated: scripts/pair.py scripts/pandoras.py
	python $< $(L2FILE) $(CONCTMPL) $(METCRO3DTMPL) $(OUTTMPL) && date > output/updated

$(FIGOUT): scripts/plot.py output/updated
	python scripts/plot.py $(FIGOUT) $(FIGIN)

$(L2FILE):
	wget -r --continue $@

clean: cleandata cleanfigs

cleandata:
	rm -f output/*

cleanfigs:
	rm -f figs/*
