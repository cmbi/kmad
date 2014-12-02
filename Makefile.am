bin_PROGRAMS	=	kman test_kman

shared_LDADD	=	$(BOOST_PROGRAM_OPTIONS_LIB) \
								$(BOOST_UNIT_TEST_FRAMEWORK_LIB)
								


kman_SOURCES	=	src/kman.cpp \
								src/fasta.cpp \
								src/features_profile.cpp \
								src/profile.cpp \
								src/residue.cpp \
								src/scoring_matrix.cpp \
								src/sequences.cpp \
								src/txtproc.cpp \
								src/vec_util.cpp \
								src/misc.cpp \
								src/substitution_matrix.cpp \
								src/msa.cpp


kman_LDADD = $(shared_LDADD)

test_kman_SOURCES = tests/test_vecutil.cpp \
										src/fasta.cpp \
										src/features_profile.cpp \
										src/profile.cpp \
										src/residue.cpp \
										src/scoring_matrix.cpp \
										src/sequences.cpp \
										src/txtproc.cpp \
										src/vec_util.cpp \
										src/misc.cpp \
										src/substitution_matrix.cpp \
										src/msa.cpp

test_kman_LDADD	=	$(shared_LDADD)


AM_CPPFLAGS	=	-std=c++11 \
							-pedantic \
							-Wall \
							-Werror \
							-Wno-reorder \
							$(BOOST_CPPFLAGS) \
							-Isrc/

if HAVE_CPPCHECK
cppcheck:
		$(CPPCHECK) --quiet --std=c++11 --enable=all \
								--suppress=missingIncludeSystem src
check:	cppcheck
endif


