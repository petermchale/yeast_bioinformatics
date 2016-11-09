from extract import getFasta
from auxFunctions import binList
from auxFunctions import tabListTuple

xMax = 1000
binLength = 100
promLength = 500
cutOffLinkerLength = 90

start_index = 0
end_index = 1
linker_index = 2
promLinkerLen_index = 5


def extractPromoters(chr, chr_string):
    """ parse saccharomyces_cerevisiae_chrxx.gff and extract promoters """

    fin = open('saccharomyces_cerevisiae_' + chr_string + '.gff')

    line = fin.readline()
    while line[0] == "#":
        line = fin.readline()

    features = []
    while line:
        seqid, source, feature_type, start, end, score, strand, phase, attributes \
            = line.split()
        if feature_type == 'CDS':

            attributes = attributes.split(';')
            initDict = [attribute.split('=') for attribute in attributes]
            attributes = dict(initDict)
            systematicGeneName = attributes['Parent']
            if 'orf_classification' in attributes:
                classification = attributes['orf_classification']
            else:
                classification = '.'
            if 'gene' in attributes:
                standardGeneName = attributes['gene']
            else:
                standardGeneName = '.'

            if strand == '+':
                promStart = int(start) - promLength
                promEnd = int(start)
            elif strand == '-':
                promStart = int(end)
                promEnd = int(end) + promLength

            promoter = chr[max(promStart, 0):promEnd]
            features += [(standardGeneName, systematicGeneName,
                          classification, promStart, promEnd, promoter)]
        line = fin.readline()

    fin.close()

    return features


def extractLinkers():
    """ extract linkers """

    fin = open('findLinkersFromRando.tab')
    line = fin.readline().split()

    features = []
    while line:
        linkerStart = int(line[start_index])
        linkerEnd = int(line[end_index])
        linker = line[linker_index]
        features += [(linkerStart, linkerEnd, linker)]
        line = fin.readline().split()

    fin.close()

    return features


def filterRandoData():

    linkerFeatures = extractLinkers()

    linkerLengthList = []
    maxLinkerLength = -1
    for linkerFeature in linkerFeatures:
        (linkerStart, linkerEnd, linker) = linkerFeature
        linkerLengthList += [len(linker)]
        if len(linker) > maxLinkerLength:
            maxLinkerLength = len(linker)
    print('max linker length =', maxLinkerLength)
    print('max length on histogram =', xMax)
    outFileName = 'allLinkerLengths.dat'
    fh = open(outFileName, 'w')
    for length in linkerLengthList:
        fh.write(str(length) + '\n')
    fh.close()

    binnedLinkerLengthList = binList(linkerLengthList, 0, xMax, binLength)
    outFileName = 'binnedLinkerLength.dat'
    fh = open(outFileName, 'w')
    for tuple_ in binnedLinkerLengthList:
        myStr = tabListTuple(tuple_)
        fh.write(myStr + '\n')
    fh.close()

    fg = open('chr03.fsa')
    (header, chr03) = getFasta(fg)
    fg.close()
    print("chr03 length =", len(chr03))

    promoterFeatures = extractPromoters(chr03, 'chr03')

    promoterLinkers = []
    for linkerFeature in linkerFeatures:
        (linkerStart, linkerEnd, linker) = linkerFeature
        if len(linker) > cutOffLinkerLength:
            for promoterFeature in promoterFeatures:
                (standardGeneName, systematicGeneName,
                 classification, promStart, promEnd, promoter) \
                    = promoterFeature
                if (promStart < linkerStart) and (linkerEnd < promEnd):
                    promoterLinker = (systematicGeneName, classification,
                                      linkerStart, linkerEnd, linker,
                                      len(linker))
                    promoterLinkers += [promoterLinker]
                    break

    print("# promoter linkers =", len(promoterLinkers))

    fileName = 'promoterLinkersFromRando.tab'
    fout = open(fileName, 'w')
    fileName = 'promoterLinkerLengths.dat'
    fout2 = open(fileName, 'w')
    for promoterLinker in promoterLinkers:
        tabbedStr = tabListTuple(promoterLinker)
        fout.write(tabbedStr + '\n')
        length = promoterLinker[promLinkerLen_index]
        fout2.write(str(length) + '\n')
    fout.close()
    fout2.close()

if __name__ == '__main__':

    filterRandoData()
