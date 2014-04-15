OBJECTS = MSA.o Sequences.o ScoringMatrix.o FeaturesProfile.o Profile.o substitutionMatrix.o misc.o txtProc.o vecUtil.o findVal.o
CLANG = clang++ -std=c++11 -stdlib=libc++
CPPFLAGS = -Os  -c  

MSA: $(OBJECTS)
	$(CLANG) -Os -o $@ $(OBJECTS) /Users/joanna/Documents/MASTER_THESIS/software/boost/stage/lib/libboost_program_options.a 
MSA.o: MSA.cpp 
	$(CLANG) $(CPPFLAGS) MSA.cpp 
Sequences.o: Sequences.h  Sequences.cpp
	$(CLANG) $(CPPFLAGS) Sequences.cpp 
ScoringMatrix.o: ScoringMatrix.h ScoringMatrix.cpp
	$(CLANG) $(CPPFLAGS) ScoringMatrix.cpp 
FeaturesProfile.o: FeaturesProfile.h FeaturesProfile.cpp
	$(CLANG) $(CPPFLAGS) FeaturesProfile.cpp
Profile.o: Profile.h Profile.cpp
	$(CLANG) $(CPPFLAGS) Profile.cpp 
substitutionMatrix.o: substitutionMatrix.h substitutionMatrix.cpp
	$(CLANG) $(CPPFLAGS) substitutionMatrix.cpp
misc.o: misc.h misc.cpp
	$(CLANG) $(CPPFLAGS) misc.cpp
txtProc.o: txtProc.h txtProc.cpp
	$(CLANG) $(CPPFLAGS) txtProc.cpp
vecUtil.o: vecUtil.h vecUtil.cpp
	$(CLANG) $(CPPFLAGS) vecUtil.cpp
findVal.o: findVal.h findVal.cpp
	$(CLANG) $(CPPFLAGS) findVal.cpp
clean: 
	rm -f *.o MSA *.h.gch
install:
	cp MSA /usr/local/bin
