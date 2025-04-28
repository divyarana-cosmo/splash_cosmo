from colossus.cosmology import cosmology
from colossus.lss import peaks

cosmology.setCosmology('planck15')
M=1e13
z=0.0
print(peaks.peakHeight(M, z))
