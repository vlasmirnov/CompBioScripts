#!/usr/bin/env python

import numpy
import sys
from argparse import ArgumentParser, SUPPRESS

'''
Thanks to Paul Zaharias for making this user-friendly
'''

def readFromFasta(filePath, removeDashes = False):
    sequences = {}
    currentTag = None

    with open(filePath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):           
                currentTag = line[1:]
                sequences[currentTag] = ""
            else :
                if(removeDashes):
                    line = line.replace("-", "")
                sequences[currentTag] = sequences[currentTag] + line        
                    
    print("Read " + str(len(sequences)) + " sequences from " + filePath + " ..")
    return sequences

def fragmentSequences(sequences, portion, fragAvg, fragStd):
    nofrags = {}
    frags = {}
   
    sequenceLengths = [len(sequences[t]) for t in sequences]
    medianLength = numpy.median(sequenceLengths)
    print("Median length {}..".format(medianLength))
    fragLength = medianLength * fragAvg
     
    taxa = list(sequences.keys())
    numpy.random.shuffle(taxa)
    numFragment = int(len(taxa)*portion)
    print("Fragmenting {} out of {}..".format(numFragment, len(taxa)))
   
    for i in range(len(taxa)):
        taxon = taxa[i]
        seq = sequences[taxon]
        if i < numFragment:            
            length = int(round(numpy.random.normal(fragLength, fragStd)))
            print("--fragment length {} out of {}..".format(length, len(seq)))
            index = numpy.random.randint(len(seq)-length+1)
            frags[taxon] = seq[index:index+length]
        else:
            nofrags[taxon] = seq

    return nofrags, frags

def saveTrueFragAligns(frags, sequences):    
    align = dict(sequences)
    for tag in frags:
        frag = str(frags[tag])
        full = str(align[tag])
        idx = full.replace("-","").find(frag)
       
        found = 0
        newString = []
        for i in range(len(full)):
            if full[i] != '-':
                if found >= idx and found < idx + len(frag):
                    newString.append(full[i])
                    assert full[i] == frag[found-idx]
                else:
                    newString.append("-")
                found = found + 1
            else:
                newString.append("-")
        align[tag] = ''.join(newString)
    return align

def write_fasta(seqs, fasta_file):
    with open(fasta_file, 'w') as f:
        for gid, gseq in seqs.items():
            f.write('>{}\n'.format(gid))
            f.write('{}\n'.format(gseq))

def main(args):
    my_true = readFromFasta(args.input_true)
    my_unaligned = {key: val.replace('-','') for key,val in my_true.items()}
    nofrags, frags = fragmentSequences(my_unaligned, args.portion_frag, args.average,args.stdev)
    align = saveTrueFragAligns(frags, my_true)
    write_fasta(align, args.output)
    if args.output_full is not None:
        write_fasta(nofrags, args.output_full)
    elif args.output_full is None:
        pass
    if args.output_frags is not None:
        write_fasta(frags, args.output_frags)
    elif args.output_frags is None:
        pass

if __name__ == "__main__":
    
    parser = ArgumentParser(description='A script to get a [partially] fragmented alignment from a given alignment.',add_help=False)

    required = parser.add_argument_group('required arguments')
    
    required.add_argument("-i", "--input_true", type=str,
                        help='input fasta of alignment.', required=True)

    required.add_argument("-p", "--portion_frag", type=float,
                        help='Proportion of sequences to be fragmented.', required=True)

    required.add_argument("-a", "--average", type=float, 
                        help='Average portion of the original median sequence length from which to draw the normal distribution.', required=True)

    required.add_argument("-s", "--stdev", type=float,
                        help='Standard deviation to the average to draw the normal distribution.', required=True)
    
    required.add_argument("-o", "--output", type=str,
                        help='Output path to alignement with fragmented sequences.', required=True)    
    
    optional = parser.add_argument_group('optional arguments')
    
    optional.add_argument("-ol", "--output_full", type=str,
                        help='Output path to fasta with unaligned full length sequences only.', required=False)

    optional.add_argument("-of", "--output_frags", type=str,
                        help='Output path to fasta with unaligned fragments.', required=False)
    
    optional.add_argument(
        '-h',
        '--help',
        action='help',
        default=SUPPRESS,
        help='show this help message and exit'
    )   

    main(parser.parse_args())
