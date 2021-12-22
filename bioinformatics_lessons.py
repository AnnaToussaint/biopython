#!/usr/bin/env python3

from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print(my_seq)

my_seq.complement()
#parsing different file formats of sequence data retrived from NCBI
from Bio import SeqIO
for seq_record in SeqIO.parse("./data/ls_orchid.fasta", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

for seq_record in SeqIO.parse("./data/ls_orchid.gbk" , "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

#difference between seq objects and strings is that different methods
#are applied to seq, though it supports some of the same methods
#seen in regular strings

my_seq = Seq("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

   # print(f"{index} {letter}")
   # an alternative way to write print statement
   #print(len(my_seq))

#access elements of the sequence, Python counts from zero index
print(my_seq[0]) #first letter in sequence
print(my_seq[2]) #third "          "
print(my_seq[-1]) #last "          "

#using .count() method
mypolyAtail = "AAAA".count("AA")
print(mypolyAtail)
mynewseq = Seq("AAAA").count("AA")
print(mynewseq)

my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
len(my_seq)

my_seq.count("G")
my_seq.count("A")

#calculating GC%
x = 100 * float(my_seq.count("G") + my_seq.count("C")) / len(my_seq)
print(x)
#alternatively, you can use Bio.SeqUtils module that has GC functions

from Bio.SeqUtils import GC
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
y = GC(my_seq)
print(y)

from Bio.Seq import Seq
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
z = my_seq[4:12]
print(z)

#get the first, second and third codon positions of a DNA sequence
my_seq[0::3]
my_seq[1::3]
my_seq[2::3]

#Turn Seq into str
str(my_seq)

#Using a placeholder for the Seq object when using Python string formatting or interpolation perator.
fasta_format_str = ">Name\n%s\n" % my_seq
print(fasta_format_string)
