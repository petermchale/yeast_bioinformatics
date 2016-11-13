from math import log  # natural log


def getFasta(fh):
    """ read a FASTA entry from the file handle fh """
    
    header = fh.readline()
    if not header:
        return header  # EOF has been found
    if header[0] != '>':
        return None  # entry is not in fasta format
    
    seq = ""
    line = fh.readline()
    while line:
        if line[0] == '>':
            # go back to the start of the header line in preparation for
            # reading next fasta entry
            fh.seek(-len(line), 1)
            break
        # remove leading and trailing numbers and white space
        line = line.strip(' 1234567890\t\n') 
        str_list = line.split()  # split into strings based on white space
        seq += "".join(str_list)  # join together all strings and add to seq
        line = fh.readline()
        
    return [header[:-1], seq]


def createEnergyMatrix(input_filename):
    """ compute the energy matrix for Gal4 using its published affinity data """

    siteLen = 17
    halfSiteLen = 8

    matrix = []
    for pos in range(siteLen):
        matrix += [{'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0}]

    fin = open(input_filename)
    lines = fin.readlines()
    count = 0
    for line in lines:
        line = line.split()
        initDictRight = {'A': float(line[0]), 'C': float(line[1]),
                         'G': float(line[2]), 'T': float(line[3])}
        initDictLeft = {'T': float(line[0]), 'G': float(line[1]),
                        'C': float(line[2]), 'A': float(line[3])}
        matrix[count+halfSiteLen] = dict(initDictRight)
        matrix[halfSiteLen-count] = dict(initDictLeft)
        count += 1
    fin.close()

    for pos in range(siteLen):
        if pos == halfSiteLen:
            for bp in 'ACGT':
                matrix[pos][bp] = -log(matrix[halfSiteLen][bp])
        else:
            for bp in 'ACGT':
                matrix[pos][bp] = -0.5*log(matrix[pos][bp])

    for pos in range(siteLen):
        minEnergy = 1000.0 
        for bp in 'ACGT':
            if matrix[pos][bp] < minEnergy:
                minEnergy = matrix[pos][bp]
        for bp in 'ACGT':
            matrix[pos][bp] -= minEnergy
            
    return matrix
