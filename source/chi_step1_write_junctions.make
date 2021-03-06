# make with make -f chi_step1_write_junctions.make

CC=g++
CFLAGS=-c -Wall -O3
OFLAGS = -Wall -O3
LDFLAGS= -Wall
SOURCES=chi_step1_write_junctions.cpp LSDMostLikelyPartitionsFinder.cpp LSDIndexRaster.cpp LSDRaster_lw.cpp LSDFlowInfo.cpp LSDChannelNetwork_lw.cpp LSDIndexChannel.cpp LSDChannel.cpp LSDIndexChannelTree.cpp LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=chi1_write_junctions.exe

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
