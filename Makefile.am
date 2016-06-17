bin_PROGRAMS	=	kmad test_kmad

shared_LDADD	=	$(BOOST_FILESYSTEM_LIB) \
			$(BOOST_PROGRAM_OPTIONS_LIB) \
			$(BOOST_REGEX_LIB) \
			$(BOOST_SYSTEM_LIB) \
			$(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
			-lconfig++

shared_SOURCES	=	src/compare_doubles.cpp \
			src/fasta.cpp \
			src/feature_scores.cpp \
			src/f_config.cpp \
			src/outfile.cpp \
			src/optimizer.cpp \
			src/profile.cpp \
			src/scoring_matrix.cpp \
			src/seq_data.cpp \
			src/msa.cpp

include_HEADERS = src/compare_doubles.h \
		  src/fasta.h \
		  src/feature_scores.h \
		  src/f_config.h \
		  src/outfile.h \
		  src/optimizer.h \
		  src/profile.h \
		  src/scoring_matrix.h \
		  src/seq_data.h \
		  src/substitution_matrix.h \
		  src/msa.h


kmad_SOURCES = $(shared_SOURCES) \
	       src/kmad.cpp

kmad_LDADD = $(shared_LDADD)

test_kmad_SOURCES = $(shared_SOURCES) \
		    tests/test_compare_doubles.cpp \
		    tests/test_config.cpp \
		    tests/test_fasta.cpp \
		    tests/test_seq_data.cpp \
		    tests/test_optimizer.cpp \
		    tests/test_profile.cpp \
		    tests/test_feature_scores.cpp \
		    tests/test_scoring_matrix.cpp \
		    tests/test_msa.cpp \
		    tests/test_outfile.cpp

EXTRA_DIST = tests/test_conf_file.cfg \
	     tests/test_conffile_pattern.cfg \
	     tests/test_conf_file_wrongformat.cfg \
	     tests/test_conf_file_settingnotfound.cfg \
	     tests/invalid_codon.fasta.7c \
	     tests/TAU_SPECI.fasta.7c \
	     tests/test_expected_feat.7c \
	     tests/test_expected_nofeat.1c \
	     tests/test_fasta_feat.7c \
	     tests/test_fasta_nofeat.1c \
	     tests/test_out_decoded.fasta_al \
	     tests/test_out_encoded.fasta_al \
	     tests/wrong_format.map \
	     tests/wrong_format2.map \
	     tests/wrong_format3.map \
	     tests/wrong_probs_format.fasta.7c \
	     tests/test.map





test_kmad_LDADD	= $(shared_LDADD)

AUTOMAKE_OPTIONS = foreign

if HAVE_CPPCHECK
cppcheck:
	$(CPPCHECK) --quiet --std=c++11 --enable=all \
		--suppress=missingIncludeSystem $(srcdir)/src
check:	cppcheck
endif
if DEBUG
AM_CXXFLAGS = -g3 -O0 \
	      -std=c++11 \
	      -pedantic \
	      -Wall \
	      -Werror \
	      -Wno-reorder \
	      $(BOOST_CPPFLAGS) \
	      -I$(srcdir)src/
else
AM_CXXFLAGS = -O3 \
	      -std=c++11 \
	      -pedantic \
	      -Wall \
	      -Werror \
	      -Wno-reorder \
	      $(BOOST_CPPFLAGS) \
	      -I$(srcdir)src/
endif
