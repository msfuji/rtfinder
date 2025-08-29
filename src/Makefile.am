bin_PROGRAMS = rttrain rtfinder 
noinst_LIBRARIES = librtfinder.a
bin_SCRIPTS = translate.rb 


librtfinder_a_SOURCES = \
genes.cpp genes.h \
fasta.cpp fasta.h \
nt.cpp nt.h \
codon.cpp codon.h \
orf.cpp orf.h


AM_CFLAGS = -g -O2 -Wall


rttrain_SOURCES = rttrain.cpp librtfinder.a 
rttrain_LDADD = librtfinder.a

rtfinder_SOURCES = rtfinder.cpp librtfinder.a
rtfinder_LDADD = librtfinder.a

