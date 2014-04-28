import sys
import os
import numpy

#input = sys.argv[1]
inputSequence = []
inputDSSP = []

aminoCode = {'A': 01,'C': 02,'D': 03,'E': 04,'F': 05,'G': 06,'H': 07,'I': 8,'K': 9,'L': 10,'M': 11,'N': 12,'P': 13,'Q': 14,'R': 15,'S': 16,'T': 17,'V': 18,'W': 19,'Y': 20}

dsspCode = {'H': 01,'G': 01,'E': 02,'B': 03,'I': 01,'S': 03,'T': 03,'C': 03,'L': 03,'_':03,'?':03}


def getInput(inputDir):
    """ Helper function for creating a list of input sequences of amino acids  """
    count = 0
    skip = False
    for root,dirs,files in os.walk(inputDir):
        for filename in files:
            for line in open(inputDir+'/'+filename,'r'):
                token = line.split(':')
                if token:
                    if token[0] not in {'RES','DSSP'}:
                        continue
                    else:
                        inp = processLine(token[1])
                        temp = list(inp)
                        if token[0] == 'RES':
                            if 'B' in inp or 'Z' in inp or 'X' in inp:
                                print "Skipping file: "+ filename
                                skip = True
                                continue
                            inputSequence.extend(inp)
                        else:
                            if skip:
                                count += 1
                                skip = False
                                continue
                            inputDSSP.extend(inp)
    print 'Skipped %d input files' % count
    return createSamples(inputSequence, inputDSSP)

def processLine(line):
    """ Helper function to split and return the input as a list """
    inp = line.strip().split(',')
    inp.remove('')
    return inp


def createSamples(res, dssp):
    """ Helper to create sliding window and corresponding y value  """
    inp = []
    labH = []
    labE = []
    labC = []
    label = []
    count = 0
    print 'Creating input profile'
    for i in range(0, len(dssp)-4):
        count += 1
#        print 'Processing sample : %d' %count
        inp.append("".join(res[i:i+5]))
        if dssp[i+2] in {'H','G','I'}:
            label.append(1.0)
            labH.append(1.0)
            labE.append(-1.0)
            labC.append(-1.0)
        else:
            labH.append(-1.0)
            if dssp[i+2] in {'E'}:            
                label.append(2.0)
                labE.append(1.0)
                labC.append(-1.0)
            else:
                label.append(3.0)
                labE.append(-1.0)
                labC.append(1.0)

    print 'Finished generating profile'
    print 'No of samples: %d' % count
    featureMatrix = createFeatureMatrix(inp)
    return {'inputs':inp,'inputUncut':inputSequence, 'label':label, 'label_H':labH, 'label_E':labE, 'label_C':labC, 'featMatrix':featureMatrix }


def createFeatureMatrix(inputSlices):
    """ Helper function to orthogonalize the input features """
    matrix = []
    
    for input in inputSlices:
        vector = []
        for char in input:
            temp = [0.0]*20
            temp[aminoCode[char]-1] = 1.0
            vector.extend(temp)
        matrix.append(vector)
    return matrix


