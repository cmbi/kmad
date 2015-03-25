# Installation
1. KMAD
KMAD can be compiled and installed from the kman directory:
  `./configure; make ; sudo make install`
2. ELM update script
To be able to use the python script to update the elm (scripts/update_elm.py)
database you need to run first:
   `pip install -r requirements.txt` 
3. convert script 
To be able to annotate predicted phosphorylations you will need to install
NetPhos. You can download it here:
http://www.cbs.dtu.dk/services/NetPhos/
Install it according to the provided instructions - make sure you can run it
simply by calling `netphos`.
The script for converting the input will run anyway without NetPhos, only you
won't have the predicted phosphorylations.


# Running

The simplest alignment of fasta sequences can be performed with the command:

  `kmad -i input.fasta -o output_prefix`

To perform alignment with features you need to first create a `.7c` file with
  `convert.py input_filename output_filename`.
The convert.py script is in the 'scripts' directory.

KMAN command for aligning with features is:
(you have to use the '-c 7' flag when running on the '.7c' input)

  `kmad -i input.7c -o output_prefix -p 10 -d 3 -m 3 -c 7`

For more detailed information try `kman --help`

To update the elm database (which you don't need to run very often)
simply execute:
  `update_elm.py`

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
pattern - regular expression by which feature X will be assigned to all
          sequences
positions - positions on which the feature is present:
            seq - sequence number (starting at 1); type = integer,
            pos - position numbers (starting at 1); type = list of integers;
            OBLIGATORY
