#!/usr/bin/env python2.7

from operator import sub
import unittest
from dnaseqlib import *

### Utility classes ###

# Maps integer keys to a set of arbitrary values.
class Multidict:
    # Initializes a new multi-value dictionary, and adds any key-value
    # 2-tuples in the iterable sequence pairs to the data structure.
    def __init__(self, pairs=[]):
        self.data = {}
        for tup in pairs:
            self.put(tup[0], tup[1])
             
    # Associates the value v with the key k.
    def put(self, k, v):
        if k in self.data:
            self.data[k].append(v)
        else:
            L = []
            self.data[k] = L
            self.data[k].append(v)
        
    # Gets any values that have been associated with the key k; or, if
    # none have been, returns an empty sequence.
    def get(self, k):
        if k in self.data:
            return self.data[k]
        else : 
            return []

# Given a sequence of nucleotides, return all k-length subsequences
# and their hashes.  (What else do you need to know about each
# subsequence?)
def subsequenceHashes(seq, k):
    subseq = ""
    for i in range(k):
        subseq += seq.next()
    print(f"subseq : {subseq}")
    arrSubSeq = list(subseq)
    rHash = RollingHash(subseq)
    # hashes = {}
    position = k
    while True:
        try:
            hash = rHash.current_hash()
            yield (hash , (subseq, position - k))
            char = seq.next()
            
            # if hash in hashes:
            #     hashes[hash].append((subseq, position - k))
            # else :
            #     L = []
            #     hashes[hash] = L
            #     hashes[hash].append((subseq, position - k))
            
            rHash.slide(arrSubSeq[0], char)
            arrSubSeq.append(char)
            del arrSubSeq[0]
            subseq = "".join(arrSubSeq)
            position += 1
        except (StopIteration):
            break
    # print(f"hashes are : {hashes}")
    # return hashes



# def subsequenceHashes(seq, k):
#     print(f"seq : {seq}")
#     subseq = seq[0 : k]
#     arrSubSeq = list(subseq)
#     rHash = RollingHash(subseq)
#     hashes = {}
#     position = k
#     for char in sequenceGenerator(seq[k:]):
#         hash = rHash.current_hash()
#         if hash in hashes:
#             hashes[hash].append((subseq , position - k))
                        
#         else:
#             L = []
#             hashes[hash] = L
#             hashes[hash].append((subseq, position - k))
            
#         rHash.slide(arrSubSeq[0], char)
#         arrSubSeq.append(char)
#         del arrSubSeq[0]
#         subseq = "".join(arrSubSeq)
#         position += 1
#     return hashes

def sequenceGenerator(seq):
    for char in seq:
        yield char

# Similar to subsequenceHashes(), but returns one k-length subsequence
# every m nucleotides.  (This will be useful when you try to use two
# whole data files.)
def intervalSubsequenceHashes(seq, k, m):
    raise Exception("Not implemented!")

# Searches for commonalities between sequences a and b by comparing
# subsequences of length k.  The sequences a and b should be iterators
# that return nucleotides.  The table is built by computing one hash
# every m nucleotides (for m >= k).
def getExactSubmatches(a, b, k, m):
    multiDictA = Multidict( subsequenceHashes(a, k))
    multiDictB = Multidict(subsequenceHashes(b, k))

    for hashA, (subseqA, positionA) in subsequenceHashes(a, k):
        if len(multiDictB.get(hashA)) != 0:
            for val in multiDictB.get(hashA):
                print(val)
                # if val[0] == subseqA:
                #     yield(positionA, val[1]) 







# def getExactSubmatches(a, b, k, m):
#     hashesA = subsequenceHashes(a, k)
#     hashesB = subsequenceHashes(b, k)
#     # print(f"hashesA : {hashesA}")
#     result = []
#     for hashKey in hashesA:
#         if hashKey in hashesB:
#             for valsA in hashesA[hashKey]:
#                 for valsB in hashesB[hashKey]:
#                     if valsA[0] == valsB[0]:
#                         result.append((valsA[1], valsB[1]))

#     return result 

if __name__ == '__main__':
    if len(sys.argv) != 4:
        # print 'Usage: {0} [file_a.fa] [file_b.fa] [output.png]'.format(sys.argv[0])
        print (f"Usage: {sys.argv[0]} [file_a.fa] [file_b.fa] [output.png]")
        sys.exit(1)

    # The arguments are, in order: 1) Your getExactSubmatches
    # function, 2) the filename to which the image should be written,
    # 3) a tuple giving the width and height of the image, 4) the
    # filename of sequence A, 5) the filename of sequence B, 6) k, the
    # subsequence size, and 7) m, the sampling interval for sequence
    # A.
    compareSequences(getExactSubmatches, sys.argv[3], (500,500), sys.argv[1], sys.argv[2], 8, 100)
