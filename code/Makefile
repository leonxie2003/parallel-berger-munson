BM_SEQ=bm_seq
BM_PAR=bm_par

COMMON_OBJS=parse_fasta.o
BM_SEQ_OBJS=bm_seq.o $(COMMON_OBJS)
BM_PAR_OBJS=bm_par.o $(COMMON_OBJS)

CXX = g++
CXXFLAGS = -Wall -O3 -std=c++17 -m64 -I.

all: $(BM_SEQ) $(BM_PAR)

$(BM_SEQ): $(BM_SEQ_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(BM_SEQ_OBJS)

$(BM_PAR): $(BM_PAR_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(BM_PAR_OBJS)

%.o: $.cpp $.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -rf *.o $(BM_SEQ) $(BM_PAR)
