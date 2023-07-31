PROGRAM         = JCorrAnalysis

LIBDIRARCH      = lib

CXX             = g++
CXXFLAGS        = -Wall -fPIC #-ggdb

LD              = g++
SOFLAGS         = -shared -Wno-deprecated
CXXFLAGS        += $(shell root-config --cflags)
LIBS            = $(shell root-config --libs)

#INCS += $(shell $(FASTJET)/bin/fastjet-config --cxxflags )
#LIBS += $(shell $(FASTJET)/bin/fastjet-config --libs --plugins ) $(INCS)

#LIBS	       += -L$(PYTHIA8)/$(LIBDIRARCH) -lpythia8 -llhapdfdummy
#INCS           += -I$(PYTHIA8)/include 
CXXFLAGS       += $(INCS)

HDRSDICT = src/JBaseTrack.h src/JBaseEventHeader.h
HDRS += $(HDRSDICT) src/JTreeDataManager.h \
	src/JHistos.h \
	src/AliAODPid.h \
	src/AliAODRedCov.h \
	src/AliAODTrack.h \
	src/AliAODVertex.h \
	src/AliBasicParticle.h \
	src/AliCFFrame.h \
	src/AliCFGridSparse.h \
	src/AliExternalTrackParam.h \
	src/AliLog.h \
	src/AliPID.h \
	src/AliTPCdEdxInfo.h \
	src/AliUEHist.h \
	src/AliVMisc.h \
	src/AliVParticle.h \
	src/AliVTrack.h \
	src/AliVVertex.h \
	src/AliCFContainer.h \
	src/AliUEHistograms.h \
	nanoDict.h

SRCS 		= $(HDRS:.h=.cxx)
OBJS		= $(HDRS:.h=.o)

%.o: %.cxx %.h
	#$(COMPILE.C) $(OUTPUT_OPTION) $(INCS) -c $<
	$(COMPILE.C) $(OUTPUT_OPTION) $<


$(PROGRAM):     $(OBJS) $(PROGRAM).C 
				$(CXX) -L$(PWD) $(PROGRAM).C $(CXXFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)
				@echo "$(PROGRAM) done"

.PHONY : clean debug

clean : 
	@echo cleaning up
	rm -f $(OBJS) core *Dict*  $(PROGRAM).o $(PROGRAM) || echo "already cleaned"

cl:  clean $(PROGRAM)

debug:
	echo $(OBJS)

nanoDict.cc: $(HDRSDICT)
		@echo "Generating dictionary ..."
		@rm -f nanoDict.cc nanoDict.hh nanoDict.h
		@rootcling nanoDict.cc -c -Done $(HDRSDICT) src/JBaseEventHeader.h
