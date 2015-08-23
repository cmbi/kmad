#!/usr/bin/python
'''
annotate a fasta sequence for KMAD with structural information

requirements:
 - SwissProt database:
     - fasta and text files split in separate entries)
     - formatted for blast
 - DSSP
 - blastp

input:
    - a regular fasta file
    or
    - a fasta.7c file (converted with IDP converter fr KMAD) - then structural
    annotations will be added to the IDP annotations

fasta headers:
    1. if you want to use information from specific PDB files then you need to
       include the PDB id in the header (where normally the uniprot ID is):
           >xxx|PDB_ID|xxx xxx xxx
       the important part is the placement of the PDB id between two pipe ('|')
       characters - the length and content of the xxx does not matter
       (unless there are no more pipes before the PDB id)
    2. If you don't know (or don't care) which PDB files you want to use,
       but you do know the Uniprot accession numbers of sequences then
       place them like the PDB_ID
       above (so just like in any fasta file from Uniprot)
       - this will speed up the process of annotation because the script doesn't
       have to run blast then. Remember to use th accession number, and not
       the entry name (e.g. for crambin use P01542 and not CRAM_CRAAB)
    3. If you didn't specify any of the aforementioned ids, then the program
       will attempt to run blast against swissprot. If the best hit has
       identity
       higher than 90%, it aligns the query sequence with it, and 'copies'
       the 2ndary structure annotation for the identical positions

'''
import argparse
import hjson
import logging
import modeller as m
import os
import re
import sys
import subprocess
import tempfile

from jsonlibconfig import encoder


SWISS_DAT = "/home/joanna/data/swissprot_dat/uniprot_dat/"
SWISS_FASTA = "/home/joanna/data/swiss_fasta/uniprot_fasta/"
SWISS_BLAST = "/home/joanna/data/uniprot_sprot"
PDB_BLAST = "/home/joanna/data/pdb_seqres.txt"
PDB_DIR = "/mnt/cmbi4/pdb/flat"
DSSP_DIR = "/home/joanna/data/dssp/"
PDBFIND = "/mnt/cmbi4/pdbfinder/PDBFIND.TXT"
SCRIPT_PATH = os.path.realpath(__file__)
KMAD = '/'.join(SCRIPT_PATH.split('/')[:-2] + ['kmad'])


log_file = "test.log"
if os.path.exists(log_file):
    os.remove(log_file)
logging.basicConfig(filename=log_file, level=logging.DEBUG)


def parse_dssp(dssp, chain_id):
    result = {"seq": "", "strct": {'G': [], 'E': [], 'B': [],
                                   'I': [], 'T': [], 'H': [],
                                   'S': [], 'C': []}}
    in_strct_section = False
    in_chain = False
    counter = 0
    for i in dssp:
        if i.startswith("  #  RESIDUE AA "):
            in_strct_section = True
        elif in_strct_section:
            if i.split()[2] == chain_id:
                in_chain = True
        if in_chain:
            if i.split()[3].islower():
                result['seq'] += 'C'
                result['strct']['C'].append(counter)
            else:
                result['seq'] += i.split()[3]

            if i.split()[4] in result['strct'].keys():
                result['strct'][i.split()[4]].append(counter)
            counter += 1
    return result


def get_strct_from_dssp(query_seq):
    strct_elements = {'G': [], 'E': [], 'B': [], 'I': [], 'T': [], 'H': [],
                      'S': [], 'C': []}
    blast_result = run_blast(query_seq, PDB_BLAST)
    if blast_result:
        chain_id = blast_result[0].split(',')[1].split('_')[1].upper()
        pdb_id = blast_result[0].split(',')[1].split('_')[0]
        dssp_path = os.path.join(DSSP_DIR, pdb_id + '.dssp')
        if os.path.exists(dssp_path):
            with open(dssp_path) as a:
                dssp = a.read().splitlines()
            dssp_data = parse_dssp(dssp, chain_id)
            equivalent_positions = align(query_seq, dssp_data['seq'])
            for i in dssp_data['strct']:
                for j in dssp_data['strct'][i]:
                    if j in equivalent_positions.keys():
                        strct_elements[i].append(equivalent_positions[j])
    return strct_elements


def run_blast(sequence, blastdb):
    tmp_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    with tmp_file as f:
        f.write(sequence)
    out_blast = tmp_file.name + '.blastp'
    args = ["blastp", "-query", tmp_file.name, "-evalue", "1e-5",
            "-num_threads", "15", "-db", blastdb,
            "-out", out_blast, '-outfmt', '10',
            "-max_target_seqs", '10']
    try:
        subprocess.call(args)
    except subprocess.CalledProcessError as e:
        print "Error: {}".format(e.output)

    if os.path.exists(out_blast):
        with open(out_blast) as a:
            output = a.read().splitlines()
        os.remove(out_blast)
    else:
        output = []
    os.remove(tmp_file.name)
    return output


# return the closest homologue sequence (if %id
# >= 90% - otherwise return empty list)
def find_closest_hit(sequence):
    blast_result = run_blast(sequence, SWISS_BLAST)
    result = {}
    # if blast_result and float(blast_result[0].split(',')[2]) >= 70:
    if blast_result:
        seq_id = blast_result[0].split(',')[1].split('|')[2]
        fasta_path = os.path.join(SWISS_FASTA, seq_id + '.fasta')
        data_path = os.path.join(SWISS_DAT, seq_id + '.dat')
        if os.path.exists(fasta_path) and os.path.exists(data_path):
            with open(fasta_path) as a:
                result['fasta'] = a.read().splitlines()
            result['data_path'] = data_path
    return result


# position - position in the aln_sequence
# returns position in the sequence without gaps
def get_real_position_al(aln_sequence, position):
    cut = re.sub('-', '', aln_sequence[:position + 1])
    return len(cut) - 1


# position - position in the sequence without gaps
# returns position in the aligned sequence
def get_real_position_seq(aln_sequence, position):
    counter = 0
    real_pos = -1
    for i in range(len(aln_sequence)):
        if aln_sequence[i] != '-':
            counter += 1
        if position == counter - 1:
            real_pos = i
            break
    return real_pos


def align(seq1, seq2):
    tmp_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    with tmp_file as f:
        f.write('\n'.join(['>1', seq1, '>2', seq2]))
    positions = {}
    args = [KMAD, '-i', tmp_file.name,
            '-o', tmp_file.name,
            '-g', '-12', '-e', '-1.2', '-n', '-1.2', '-c', '1', '--gapped']
    try:
        subprocess.call(args)
        if os.path.exists(tmp_file.name + '_al'):
            with open(tmp_file.name + '_al') as a:
                aln = unwrap(a.read().splitlines())
            aln_seq1 = aln[1]
            aln_seq2 = aln[3]
            for i in range(len(aln_seq1)):
                if aln_seq1[i] == aln_seq2[i]:
                    real_pos2 = get_real_position_al(aln_seq2, i)
                    real_pos1 = get_real_position_seq(seq1,
                                                      get_real_position_al(
                                                          aln_seq1, i))
                    positions[real_pos2] = real_pos1
            os.remove(tmp_file.name + '_al')

    except subprocess.CalledProcessError as e:
        print "Error: {}".format(e.output)
    os.remove(tmp_file.name)
    return positions


# get 2ndary strct info based on the information from a very close homologue
# (transfer only data from identical positions)
def transfer_data_from_homologue(sequence, closest_sp):
    seq2 = unwrap(closest_sp['fasta'])[1]
    equivalent_positions = align(sequence, seq2)
    homologue_strct_elems = get_strct_from_sp(closest_sp['data_path'])
    query_strct_elems = {'G': [], 'E': [], 'B': [], 'I': [], 'T': [], 'H': [],
                         'S': [], 'C': []}
    for i in homologue_strct_elems:
        for j in homologue_strct_elems[i]:
            if j in equivalent_positions.keys():
                query_strct_elems[i].append(equivalent_positions[j])
    return query_strct_elems


def change_char(mystring, position, new_char):
    if position < len(mystring) - 1:
        new_string = mystring[:position] + new_char + mystring[position + 1:]
    else:
        new_string = mystring[:position] + new_char
    return new_string


# based on the dictionary strct_data (with lists of positions of certain 2ndary
# strct elements) encode a fasta file (only for a plain fasta file, not for 7c)
def encode(fasta, strct_data):
    encoded = []
    for i, lineI in enumerate(fasta):
        if lineI.startswith('>'):
            encoded.append(lineI)
        else:
            got_data = any(strct_data[i / 2].values())
            if got_data:
                newline = ""
                for j in range(len(lineI)):
                    codon = lineI[j] + 'AAAAAA'
                    for k in strct_data[i / 2].keys():
                        if k != 'C' and j in strct_data[i / 2][k]:
                            codon = change_char(codon, 1, k)
                    if j in strct_data[i / 2]['C']:
                        codon = change_char(codon, 4, 's')
                    newline += codon
                encoded.append(newline)
            else:
                encoded.append(''.join([j + 'AAAAAA' for j in lineI]))
    return encoded


# get 2ndary strct annotations from swiss prot
def get_strct_from_sp(sp_path):
    with open(sp_path) as a:
        dat_file = a.read().splitlines()
        strct_elements = {'G': [], 'E': [], 'B': [], 'I': [], 'T': [], 'H': [],
                          'S': [], 'C': []}
        strct_dict = {'HELIX': 'H', 'TURN': 'T', 'STRAND': 'S',
                      'DISULFID': 'C'}
    for i in dat_file:
        if (i.startswith('FT') and len(i.split()) > 1
                and i.split()[1] in strct_dict.keys()):
            if i.split()[1] != 'DISULFID':
                start = int(i.split()[2]) - 1
                end = int(i.split()[3])
                strct_code = strct_dict[i.split()[1]]
                strct_elements[strct_code].extend(range(start, end))
            else:
                cys1 = int(i.split()[2]) - 1
                cys2 = int(i.split()[3]) - 1
                strct_elements['C'].extend([cys1, cys2])
    return strct_elements


def unwrap(alignment):
    new = []
    for i in alignment:
        if i.startswith('>'):
            new.append(i)
            new.append("")
        else:
            new[-1] += i
    return new


def decode(fasta_7c):
    fasta = []
    for i in fasta_7c:
        if i.startswith('>'):
            fasta.append(i)
        elif not i.startswith('#'):
            fasta.append(''.join(i[::7]))
        else:
            break
    return fasta


def get_seq_from_pdb(pdb_id, chain_id):
    with open(PDBFIND) as a:
        pdbfind = a.read().splitlines()
    in_prot = False
    in_chain = False
    sequence = ""
    for i in pdbfind:
        if i.startswith('ID'):
            if not in_prot and i.split()[-1] == pdb_id.upper():
                in_prot = True
            elif in_prot:
                break
        elif in_prot and i.startswith("Chain") and i.split()[-1] == chain_id:
            in_chain = True
        elif in_chain and i.startswith(" Sequence"):
            sequence = i.split()[-1]
            break
    return sequence


def merge_dicts(dict1, dict2):
    dict_tmp = dict1.copy()
    dict_tmp.update(dict2)
    return dict_tmp


def annotate_secondary_structure(fasta, output_name):
    strct_data = []
    for i in range(0, len(fasta), 2):
        strct_elements = {'G': [], 'E': [], 'B': [], 'I': [], 'T': [], 'H': [],
                          'S': [], 'C': []}
        closest_sp = find_closest_hit(fasta[i + 1])
        if closest_sp:
            strct_elements = transfer_data_from_homologue(fasta[i + 1],
                                                          closest_sp)
        print strct_elements
        strct_elements = merge_dicts(strct_elements,
                                     get_strct_from_dssp(fasta[i + 1]))
        print strct_elements
        strct_data.append(strct_elements)
    out_fasta = encode(fasta, strct_data)

    out = open(output_name, 'w')
    out.write('\n'.join(out_fasta))
    out.close()


def get_structure_data(query_seq, pdb_fastas):
    result = {'id': '',
              'chain_id': '',
              'pdb_path': '',
              'pdb_seq': '',
              'positions_map': {}
              }
    blast_result = run_blast(query_seq, PDB_BLAST)
    if blast_result:
        result['chain_id'] = blast_result[0].split(',')[1].split('_')[1]
        result['id'] = blast_result[0].split(',')[1].split('_')[0]
        result['pdb_seq'] = get_seq_from_pdb(result['id'], result['chain_id'])
        result['pdb_path'] = os.path.join(PDB_DIR,
                                          'pdb' + result['id'].lower() + '.ent')
        result['positions_map'] = align(result['pdb_seq'], query_seq)
    return result


# run structure alignment and return a dictionary of equivalent positions
def mk_strct_al_modeller(strct_data1, strct_data2):
    _stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w')

    tmp_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    env = m.environ()

    aln = m.alignment(env)
    code1 = 'pdb' + strct_data1['id']
    code2 = 'pdb' + strct_data2['id']
    chain1 = strct_data1['chain_id']
    chain2 = strct_data2['chain_id']
    env.io.atom_files_directory = ['.', PDB_DIR]
    result = {}
    try:
        for (code, chain) in ((code1, chain1), (code2, chain2)):
            mdl = m.model(env, file=code, model_segment=('FIRST:'+chain,
                                                         'LAST:'+chain))
            aln.append_model(mdl, atom_files=code, align_codes=code+chain)

        for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False,
                                             True),
                                            ((1., 0.5, 1., 1., 1., 0.), False,
                                             True),
                                            ((1., 1., 1., 1., 1., 0.), True,
                                             False)):
            r = aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                           rr_file='$(LIB)/as1.sim.mat', overhang=30,
                           gap_penalties_1d=(-450, -50),
                           gap_penalties_3d=(0, 3), gap_gap_score=0,
                           gap_residue_score=0,
                           alignment_type='tree', # If 'progresive', the tree is not
                                                  # computed and all structures will be
                                                  # aligned sequentially to the first
                           #ext_tree_file='1is3A_exmat.mtx', # Tree building can be avoided
                                                             # if the tree is input
                           feature_weights=weights, # For a multiple sequence alignment only
                                                    # the first feature needs to be non-zero
                           improve_alignment=True, fit=True, write_fit=False,
                           write_whole_pdb=whole, output='ALIGNMENT QUALITY')
        if r.qscorepct > 70:
            aln.write(file=tmp_file.name, alignment_format='FASTA')
            with open(tmp_file.name) as a:
                alignment = unwrap(a.read().splitlines())

            for i in range(len(alignment[1])):
                if alignment[1] != '-' and alignment[3] != '-':
                    pos1 = get_real_position_al(alignment[1], i)
                    pos2 = get_real_position_al(alignment[3], i)
                    result[pos1] = pos2
    except:
        print 'Modeller failed'
    sys.stdout.close()
    sys.stdout = _stdout
    return result


def mk_strct_alignment(strct_data1, strct_data2):
    args = ['TMalign', strct_data1['pdb_path'], strct_data2['pdb_path']]
    result = {}
    output = []
    try:
        output = subprocess.check_output(args).splitlines()
    except subprocess.CalledProcessError as e:
        print "Error: {}".format(e.output)

    if output:
        seq1 = output[-4]
        seq2 = output[-2]
        al = output[-3]
        # if (al.find(':' * (len(al) / 4)) != -1
        if len(seq1) == len(seq2) and len(seq1) == len(al):
            for i in range(len(seq1)):
                if al[i] == ':':
                    pos1 = get_real_position_al(seq1, i)
                    pos2 = get_real_position_al(seq2, i)
                    result[pos1] = pos2
        else:
            print "bad alignment {}".format(al), al.count(':'), len(al)
    return result


def reverse_dict(some_dict):
    new_dict = {}
    for i in some_dict:
        new_dict[some_dict[i]] = i
    return new_dict


# return equivalent positions between seq and seq2
def find_equivalent_positions(seq1_str1, seq2_str2, str1_str2):
    str1_seq1 = reverse_dict(seq1_str1)
    str2_seq2 = reverse_dict(seq2_str2)

    seq1_seq2 = {}

    for i in str1_str2:
        i_str1 = i
        i_str2 = str1_str2[i]
        i_seq1 = ''
        i_seq2 = ''
        if i_str1 in str1_seq1.keys():
            i_seq1 = str1_seq1[i_str1]
        if i_str2 in str2_seq2.keys():
            i_seq2 = str2_seq2[i_str2]
        if i_seq1 and i_seq2:
            seq1_seq2[i_seq1] = i_seq2
    return seq1_seq2


def make_conf_dict(query_seq, eq_positions, al_score):
    settings_dict = {"feature_settings": {"usr_features": []}}
    all_eqs = []
    for i in eq_positions:
        all_eqs.extend(eq_positions[i].keys())
    for pos in range(len(query_seq)):
        curr_feature_positions = []
        if pos in all_eqs:
            for seq_no in eq_positions:
                if pos in eq_positions[seq_no].keys():
                    curr_feature_positions.append(
                        {"seq": seq_no + 1,
                         "pos": tuple([eq_positions[seq_no][pos] + 1])})
        if curr_feature_positions:
            curr_feature_positions.append({"seq": 1, "pos": tuple([pos + 1])})
            single_feat = {"name": "pos{}".format(pos + 1),
                           "add_score": al_score,
                           "subtract_score": 0,
                           "add_features": tuple(["pos{}".format(pos + 1)]),
                           "add_tags": [],
                           "add_exceptions": [],
                           "subtract_features": [],
                           "subtract_tags": [],
                           "subtract_exceptions": [],
                           "subtract_features": [],
                           "pattern": '',
                           "positions": tuple(curr_feature_positions)}
            settings_dict["feature_settings"]["usr_features"].append(
                single_feat)
    return settings_dict


def dict_to_cfg(data_dict):
    data_h = hjson.dumps(data_dict)
    text = re.sub('\}', '};', re.sub('\]', ');', re.sub('\[', '(', data_h)))
    text = text[:-2].lstrip('{')
    return text


def write_conf_file(query_seq, eq_positions, output_conf, al_score):
    data_dict = make_conf_dict(query_seq, eq_positions, al_score)
    indent = 2
    outtxt = encoder.dumps(data_dict, indent=indent)
    out = open(output_conf, 'w')
    out.write(outtxt)
    out.close()


def structure_alignment_conf(fasta, output_conf, al_score):
    with open(PDB_BLAST) as a:
        pdb_fastas = unwrap(a.read().splitlines())
    query_sequence = fasta[1]
    query_strct_data = get_structure_data(query_sequence, pdb_fastas)
    eq_positions = {}
    for i in range(2, len(fasta)):
        if not fasta[i].startswith('>'):
            strct_data = get_structure_data(fasta[i], pdb_fastas)
            eq_seq1_seq2 = {}
            if strct_data['id']:
                strct_al = mk_strct_al_modeller(query_strct_data, strct_data)
                eq_seq1_seq2 = find_equivalent_positions(
                    query_strct_data['positions_map'],
                    strct_data['positions_map'],
                    strct_al)
            eq_positions[i / 2] = eq_seq1_seq2
    write_conf_file(query_sequence, eq_positions, output_conf, al_score)


def annotate(fasta_in, output_name, output_conf, al_score):
    with open(fasta_in) as a:
        fasta = unwrap(a.read().splitlines())
    annotate_secondary_structure(fasta, output_name)
    structure_alignment_conf(fasta, output_conf, al_score)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate fasta sequences with'
                                                 ' structural information')
    parser.add_argument('fasta_in')
    parser.add_argument('output_filename')
    parser.add_argument('output_conf')
    parser.add_argument('--al_score', type=int)
    args = parser.parse_args()
    annotate(args.fasta_in, args.output_filename, args.output_conf,
             args.al_score)
