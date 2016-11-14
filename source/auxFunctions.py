from math import floor
import sys
import pandas as pd


def calcSiteEnergyWithMatrix(seq, matrix):
    """ calculate the energy of seq using an energy matrix """

    if len(seq) != len(matrix):
        sys.stderr.write('seq and matrix not the same size!\n')
        sys.exit(-1)

    energy = 0.0 
    for pos in range(len(seq)):
        energy += matrix[pos][seq[pos]]

    return energy


def calcEnergyListWithMatrix(dna, matrix):
    """ calculate an energy for each site in DNA """

    siteLen = len(matrix)
    if siteLen > len(dna):
        sys.stderr.write('site longer than dna to be scanned\n')
        sys.exit(-1)
        
    energyList = []
    for start in range(len(dna)-siteLen+1):
        end = start + siteLen
        seq = dna[start:end]
        seqEnergy = calcSiteEnergyWithMatrix(seq, matrix)
        energyList += [(seqEnergy, start, end)]

    return pd.DataFrame(data=energyList, columns=['TF-DNA binding energy', 'binding-site start position', 'binding-site end position'])


def calc_bin(num, binWidth):
    """ determine which bin num belongs to """

    bin_ = floor(num/binWidth)
    bin_ = int(bin_)
    return bin_


def binList(xList, xMin, xMax, binWidth):
    """ bin list of energies, probabilities, etc """

    numberBins = calc_bin(xMax-xMin, binWidth) + 1

    binnedxList = []
    for i in range(numberBins):
        xBin = xMin + (0.5+i)*binWidth
        binnedxList += [(xBin, 0)]

    for x in xList:
        binNum = calc_bin(x-xMin, binWidth)
        xBin, xFreq = binnedxList[binNum]
        binnedxList[binNum] = xBin, xFreq+1

    return binnedxList


def tabListTuple(mixedListTuple):
    """ create tab-delimited string from a list/tuple """

    return "\t".join([str(item) for item in mixedListTuple]) 
