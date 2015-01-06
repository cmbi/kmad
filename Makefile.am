bin_PROGRAMS	=	kman test_kman

shared_LDADD	=	$(BOOST_FILESYSTEM_LIB) \
								$(BOOST_PROGRAM_OPTIONS_LIB) \
								$(BOOST_REGEX_LIB) \
								$(BOOST_SYSTEM_LIB) \
								$(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
								-lconfig++

shared_SOURCES	=	src/fasta.cpp \
									src/features_profile.cpp \
									src/f_config.cpp \
									src/profile.cpp \
									src/scoring_matrix.cpp \
									src/txtproc.cpp \
									src/vec_util.cpp \
									src/substitution_matrix.cpp \
									src/msa.cpp


kman_SOURCES	=	$(shared_SOURCES) \
								src/kman.cpp


kman_LDADD = $(shared_LDADD)

test_kman_SOURCES = $(shared_SOURCES) \
										tests/test_config.cpp \
										tests/test_fasta.cpp \
										tests/test_features_profile.cpp \
										tests/test_profile.cpp

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


