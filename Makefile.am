bin_PROGRAMS	=	kman test_kman

shared_LDADD	=	$(BOOST_FILESYSTEM_LIB) \
								$(BOOST_PROGRAM_OPTIONS_LIB) \
								$(BOOST_REGEX_LIB) \
								$(BOOST_SYSTEM_LIB) \
								$(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
								-lconfig++ \
								-lgcov

shared_SOURCES	=	src/fasta.cpp \
									src/feature_scores.cpp \
									src/f_config.cpp \
									src/outfile.cpp \
									src/profile.cpp \
									src/scoring_matrix.cpp \
									src/seq_data.cpp \
									src/msa.cpp


kman_SOURCES	=	$(shared_SOURCES) \
								src/kman.cpp


kman_LDADD = $(shared_LDADD)

test_kman_SOURCES = $(shared_SOURCES) \
										tests/test_config.cpp \
										tests/test_fasta.cpp \
										tests/test_msa.cpp \
										tests/test_scoring_matrix.cpp \
										tests/test_seq_data.cpp \
										tests/test_profile.cpp \
										tests/test_feature_scores.cpp

									

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


