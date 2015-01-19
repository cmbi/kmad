make clean
./configure CXXFLAGS='--coverage'
./test_kman
gcov -r -o src/msa src/msa.cpp
lcov --capture --directory src/ --output-file coverage.info
genhtml coverage.info --output-directory out
find . -name "*gcno" | xargs -r rm
find . -name "*gcda" | xargs -r rm
find . -name "*gcov" | xargs -r rm
