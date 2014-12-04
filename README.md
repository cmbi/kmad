# Installation

KMAN can be compiled and installed from the kman directory:
./configure; make ; sudo make install

# Running

The simplest alignment of fasta sequences can be performed with the command:

    kman -i input.fasta -o output_prefix -g

To perform alignment with features you need to first create a `.7c` file with
`convertFastaTo7c.py`.

KMAN command for aligning with features is:

    kman -i input.7c -o output_prefix -p 10 -d 3 -m 3 -c 7

For more detailed information try `kman --help`


# Configuration file

You can define your own features to be included in the alignment.
This can be done with a configuration file (--conf flag). An example of 
a configuration file is provided in the KMAN directory ("conf_file_example").

Description

All of the possible fields are present in the example file, but only some of
them are opbligatory. Here's a detailed decsription of the field meanings:

name - a unique feature name; type = string; OBLIGATORY

tag - features can be assigned various categories, this is the category's name;
      type = string; optional 

add_score - score for aligning features that add points to the 
            alignment score;
            type = positive double; optional

subtract_score - score for aligning features that subtract points from the
                 alignment score; type = positive double; optional
 
add_features - names of features that add points to the alignment; 
               type = list of strings; OBLIGATORY (it can be an empty list)

add_tags - names of feature categories that add points to the alignment;
          type = list of strings; OBLIGATORY (it can be an empty list)

add_exceptions - names of exceptions to the add_tags field;
                 type = list of strings; OBLIGATORY (it can be an empty list)

subtract_features - names of features that subtract points from the alignment; 
                    type = list of strings; OBLIGATORY (it can be an empty 
                    list)

subtract_tags - names of feature categories that subtract from the alignment;
                type = list of strings; OBLIGATORY (it can be an empty list)

subtract_exceptions - names of exceptions to the subtract_tags field;
                      type = list of strings; OBLIGATORY (it can be an empty 
                      list)

positions - positions on which the feature is present:
            seq - sequence number (starting at 1); type = integer,
            pos - position numbers (starting at 1); type = list of integers;
            OBLIGATORY
