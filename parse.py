import sys
import inputParser


def run(inp):
    for line in open(inp,'r'):
                token = line.split(':')
                if token:
                    if token[0] not in {'RES','DSSP'}:
                        continue
                    else:
                        inp = inputParser.processLine(token[1])
                        temp = list(inp)
                        if token[0] == 'RES':
                            if 'B' in inp or 'Z' in inp or 'X' in inp:
                                print "Skipping file: "+ inp
                                skip = True
                                continue
                            return processLine(token[1])


def processLine(line):
    input = line.strip().replace(',','')
    print "Input: " + input
    inp = line.strip().split(',')
    inp.remove('')
    return inp

