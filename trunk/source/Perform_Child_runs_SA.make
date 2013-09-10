# make with make -f Perform_Child_runs_SA.make

CC=g++
CFLAGS=-c -Wall -O3 
OFLAGS = -Wall -O3 
LDFLAGS= -Wall
SOURCES=Perform_chi_for_Child_runs_sensitivity.cpp LSDMostLikelyPartitionsFinder.cpp LSDChiNetwork.cpp LSDStatsTools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Perform_Child_SA.exe

all: $(SOURCES) $(EXECUTABLE)

#$(EXECUTABLE): $(OBJECTS)
#	$(CC) $(OFLAGS) $(OBJECTS) $(LIBS) -o $@

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
