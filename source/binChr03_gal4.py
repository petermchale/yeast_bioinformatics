from extract import createEnergyMatrix, getFasta
from auxFunctions import calcEnergyListWithMatrix, binList, tabListTuple

E_MIN = -5.0
E_MAX = 50.0
ENERGY_BIN_WIDTH = 0.25

matrix = createEnergyMatrix()
binSizeStr = str(ENERGY_BIN_WIDTH)

print("evaluating energy list ...")
fh = open('chr03.fsa')
[header, sequence] = getFasta(fh)
energyList = calcEnergyListWithMatrix(sequence, matrix)
fh.close()
print("...energy list evaluated")
print("binning energies...")
binnedEnergyList = binList(energyList, E_MIN, E_MAX, ENERGY_BIN_WIDTH)
print("...energies binned")
outFileName = 'binChr03_Gal4_' + binSizeStr + '_.dat'
fh = open(outFileName, 'w')
for tuple_ in binnedEnergyList:
    myStr = tabListTuple(tuple_)
    fh.write(myStr+'\n')
fh.close()
