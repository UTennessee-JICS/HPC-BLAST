CC=icc
CFLAGS=-O2
CXX=icpc
CXXFLAGS=-O2 -std=c++11

default: all
all: database_sample distribute_queries get_query_stats stitch_blast_output strip_header compare_results

database_sample:
	$(CC) -msse2 -DHAVE_SSE2 -DMEXP=607 -o $@ SFMT.c $@.c -lm
distribute_queries:
	$(CXX) $(CXXFLAGS) -qopenmp -o $@ $@.cpp
get_query_stats:
	$(CC) $(CFLAGS) -o $@ $@.c
stitch_blast_output:
	$(CXX) $(CXXFLAGS) -o $@ $@.cpp
strip_header:
	$(CC) $(CFLAGS) -o $@ $@.c
compare_results:
	$(CC) $(CFLAGS) -o $@ $@.c

clean:
	rm -rf database_sample distribute_queries get_query_stats stitch_blast_output strip_header compare_results *.o
