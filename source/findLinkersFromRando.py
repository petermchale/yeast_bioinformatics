# Note: the start and stop coords in Yuan et al are not consistent with the 50mers
# they use: the coords are correct only to within about 4bps

from auxFunctions import tabListTuple
from extract import getFasta

probe_index = 0
probe_start_index = 2
probe_end_index = 3
nucl_state_index = 4

LINE_LEN = 50

fg = open('chr03.fsa')
[header, chr03] = getFasta(fg)

fin = open('concatenateRandoData_bare.tab')
fout = open('findLinkersFromRando.tab', 'w')

line = fin.readline().split()
lastProbe = int(line[probe_index])
lastNuclState = int(line[nucl_state_index])
if lastNuclState == 0:
    start = int(line[probe_start_index])
    end = int(line[probe_end_index])

line = fin.readline().split()
while line:
    currentProbe = int(line[probe_index])
    currentNuclState = int(line[nucl_state_index])

    if (lastNuclState > 0) and (currentNuclState == 0):
        start = int(line[probe_start_index])
        end = int(line[probe_end_index])
    elif (lastNuclState == 0) and (currentNuclState == 0):
        if (currentProbe-lastProbe) == 2:
            end = int(line[probe_end_index])
        else:
            sequence = chr03[start-1:end]
            length = len(sequence)
            outList = [start, end, sequence, length]
            tabbedStr = tabListTuple(outList)
            fout.write(tabbedStr+'\n')
            start = int(line[probe_start_index])
            end = int(line[probe_end_index])
    elif (lastNuclState == 0) and (currentNuclState > 0):
        sequence = chr03[start-1:end]
        length = len(sequence)
        outList = [start, end, sequence, length]
        tabbedStr = tabListTuple(outList)
        fout.write(tabbedStr+'\n')

    line = fin.readline().split()
    lastProbe = currentProbe
    lastNuclState = currentNuclState

if lastNuclState == 0:
    sequence = chr03[start-1:end]
    length = len(sequence)
    outList = [start, end, sequence, length]
    tabbedStr = tabListTuple(outList)
    fout.write(tabbedStr+'\n')

fin.close()
fout.close()

###############################################################
