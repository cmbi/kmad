OBJECTS = MSA.o Sequences.o ScoringMatrix.o SubstitutionMatrix.o Profile.o UsefulStuff.o


MSA: $(OBJECTS)
	g++ -o $@ $(OBJECTS)

MSA.o: MSA.cpp 
	g++ -c MSA.cpp 

Sequences.o: Sequences.h  Sequences.cpp
	g++ -c Sequences.cpp 

ScoringMatrix.o: ScoringMatrix.h ScoringMatrix.cpp
	g++ -c ScoringMatrix.cpp 

SubstitutionMatrix.o: SubstitutionMatrix.h SubstitutionMatrix.cpp
	g++ -c SubstitutionMatrix.cpp
Profile.o: Profile.h Profile.cpp
	g++ -c Profile.cpp 
UsefulStuff.o: UsefulStuff.h UsefulStuff.cpp
	g++ -c UsefulStuff.cpp

clean: 
	rm -f *.o MSA *.h.gch
