from thinkdsp import Wave,decorate

from DPI import hl2adb as hl2a
from DPI import changeDriver
# 0=driver mode   0=hl2a,1=local, 2=exl50,3=east,4=hl2m
changeDriver(0)
shotNum = 36964
ts0 = 150 / 1e3
te0 = 900 /1e3

ts, y, U = hl2a(shotNum, 'Mpol_14', ts0, te0,1)