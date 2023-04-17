#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main

#include "src/fasta.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test.hpp>

using namespace fasta;

BOOST_AUTO_TEST_SUITE(test_fasta)


BOOST_AUTO_TEST_CASE(test_parse_fasta)
{
    FastaData fd = parse_fasta("tests/TAU_SPECI.fasta.7c", 7);
    BOOST_CHECK_EQUAL(fd.sequences.size(), 19);
    BOOST_CHECK_EQUAL(
        sequence_to_string(fd.sequences[0]),
        "MAEPRQEFDTAEDHAEGYALLQDQEGEHGLKASPLQTPADDGPEEPVSETSDAKSTPTAEDVT"
        "APLVDERTPGEQAATQPPTDIPEGTTAEEAGIGDTPNMEDQAAGHVTQARMVSKGKEGTGSED"
        "RKAKGADSKTGTKIATPRGTAPPGQKGTANATRIPAKTTPSPKTPPGTGEPAKSGDRSGYSSP"
        "GSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSTKSRLQTAPVPMPDLKNVRSKIGST"
        "ENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLG"
        "NIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHG"
        "AEIVYKSPVVSGDTSPRHLSNVSSTGSINMVDSPQLATLADEVSASLAKQGL"
    );
    BOOST_CHECK_EQUAL(
        sequence_to_string(fd.sequences[1]),
        "MAEPRQEFDVMEDHAQGDYTLQDHEGDMEPGLKESPLQTPADDGSEEPGSETSDAKSTPTAEA"
        "EEAGIGDTSNLEDQAAGHVTQARMVSKGKDGTGPDDKKAKGADGKPGTKIATPRGAAPPGQKG"
        "QANATRIPAKTTPTPKTSPGTGESGKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVA"
        "VVRTPPKSPSAAKSRLQAAPGPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSK"
        "CGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSK"
        "IGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGS"
        "IDMVDSPQLATLADEVSASLAKQGL"
    );
    BOOST_CHECK_EQUAL(fd.probabilities.size(), 97);
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_nonexistent)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/nonexistent.fasta.7c", 7), std::invalid_argument
    );
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_invalid_probabilities_format)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/wrong_probs_format.fasta.7c", 7), std::runtime_error
    );
}

BOOST_AUTO_TEST_CASE(test_parse_fasta_invalid_codons)
{
    BOOST_CHECK_THROW(
        parse_fasta("tests/invalid_codon.fasta.7c", 7), std::runtime_error
    );
}

BOOST_AUTO_TEST_CASE(test_make_sequence)
{
  Sequence s = make_sequence("desc", "AaaaaaBbbbbbCccccc", 6);
  BOOST_CHECK_EQUAL(s.residues.size(), 3);
}

BOOST_AUTO_TEST_CASE(test_sequence_to_string)
{
  Sequence s = make_sequence("desc", "LSKAL", 1);
  std::string string_seq = sequence_to_string(s);
  std::string expected = "LSKAL";
  BOOST_CHECK_EQUAL(string_seq, expected);
}


BOOST_AUTO_TEST_CASE(test_remove_gaps) {
  fasta::Sequence s1 = fasta::make_sequence("", "ABAD-C", 1);
  fasta::Sequence s2 = fasta::make_sequence("", "ATSSAC", 1);
  fasta::Sequence s3 = fasta::make_sequence("", "--BAA-", 1);
  fasta::SequenceList s = {s1, s2 , s3};
  s = fasta::remove_gaps(s);
  std::vector<std::string> result;
  for (auto& seq : s) {
    result.push_back(fasta::sequence_to_string(seq));
  }
  std::vector<std::string> expected = {"ABADC", "ATSSAC", "BAA"};
  BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
                                expected.begin(), expected.end());
  s1 = fasta::make_sequence("", "ABADC", 1);
  s2 = fasta::make_sequence("", "ATSSAC", 1);
  s3 = fasta::make_sequence("", "BAA", 1);
  s = {s1, s2, s3};
  s = fasta::remove_gaps(s);
  result.clear();
  for (auto& seq : s) {
    result.push_back(fasta::sequence_to_string(seq));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(result.begin(), result.end(),
                                expected.begin(), expected.end());
}



BOOST_AUTO_TEST_SUITE_END()
