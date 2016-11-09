from auxFunctions import tabListTuple
import sys

probe_index = 0
chromo_index = 1
v5_index = 3

fin_probe = open('randoProbeInfo.in')
fin_calls = open('randoHMMhandCalls.in')
fout_bare = open('concatenateRandoData_bare.tab', 'w')
fout_sequence = open('concatenateRandoData_seq.tab', 'w')

# header lines
line_probe = fin_probe.readline().split()
line_calls = fin_calls.readline().split()

# first line
line_probe = fin_probe.readline().split()
line_calls = fin_calls.readline().split()

# process files
while line_probe and line_calls:
    if line_probe[probe_index] != line_calls[probe_index]:
        sys.stderr.write('probes are not the same!')
        sys.exit(-1)
    if (line_probe[chromo_index] == '3') and (line_calls[v5_index] != 'NaN'):
        newLine = line_probe[:-1] + [line_calls[v5_index]]
        tabbedStr = tabListTuple(newLine)
        fout_bare.write(tabbedStr+'\n')
        newLine = line_probe + [line_calls[v5_index]]
        tabbedStr = tabListTuple(newLine)
        fout_sequence.write(tabbedStr+'\n')
    line_probe = fin_probe.readline().split()
    line_calls = fin_calls.readline().split()

fin_probe.close()
fin_calls.close()
fout_bare.close()
fout_sequence.close()
