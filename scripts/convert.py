#!/usr/bin/python
import argparse
import logging
import httplib
import math
import os
import re
import string
import subprocess
import sys
import tempfile
import time
import urllib
import urllib2

from operator import itemgetter
import xml.etree.ElementTree as ET


def create_numbering():
    alphabet = list(string.ascii_lowercase) \
        + list(string.ascii_uppercase) \
        + [str(i) for i in range(0, 10)]
    alphabet.remove('A')
    numbering = []
    for i in alphabet:
        for j in alphabet:
            numbering += [i+j]
    return numbering


def read_fasta(fastafilename):
    with open(fastafilename) as a:
        fastafile = a.read()
    fastafile_list = fastafile.splitlines()
    result = []
    for i, elementI in enumerate(fastafile_list):
        if elementI.startswith('>') or fastafile_list[i-1].startswith('>'):
            result += [elementI]
        else:
            result[-1] += elementI
    return result


def process_slims_all_classes(classes):
    result = dict()
    for i in classes:
        lineI = i.split()
        slim_id = lineI[0]
        slim_regex = lineI[1]
        slim_probability = float(lineI[2])
        go_terms = lineI[3:]
        result[slim_id] = {"prob": slim_probability, "GO": go_terms,
                           "regex": slim_regex}
    return result


def get_id(sequence):
    if len(sequence.split('|')) >= 3:
        return sequence.split('|')[2].split(' ')[0]
    else:
        return sequence.lstrip('>')


# PFAM
def run_pfam_scan(filename):
    domain_coords = []
    domain_accessions = []
    with open(filename) as a:
        fastafile = re.sub('-', '', a.read())
    values = {'seq': fastafile, 'output': 'xml'}
    data = urllib.urlencode(values)
    try:
        request = urllib2.Request('http://pfam.xfam.org/search/sequence', data)
        reply = urllib2.urlopen(request).read()
        pfam_server_error = False
        try:
            tree = ET.fromstring(reply)
        except ET.ParseError:
            pfam_server_error = True
        if not pfam_server_error:
            result_url = tree[0][1].text
            count = 0
            finished = False
            while count < 20 and not finished:
                time.sleep(6)
                try:
                    request = urllib2.Request(result_url)
                    result = urllib2.urlopen(request).read()
                    root = ET.fromstring(result)
                    if len(root[0]):
                        for child in root[0][0][0][0]:
                            for g in child:
                                domain_accessions += [child.attrib['accession']]
                                domain_coords += [[int(g.attrib['start']),
                                                   int(g.attrib['end'])]]
                    finished = True
                except ET.ParseError:
                    pass
                count += 1
    except urllib2.HTTPError:
        print "pfam scan HTTPError"
    return [domain_coords, domain_accessions]


# NETPHOS
def run_netphos(filename):
    phosphorylations = set([])
    try:
        args = ['netphos', filename]
        netPhosOut = subprocess.check_output(args).splitlines()

        for lineI in netPhosOut:
            if len(lineI.split()) > 0 and lineI.split()[-1] == 'YES':
                phosphorylations.add(int(lineI.split()[2]))
    except subprocess.CalledProcessError:
        print "Not running netphos"
    return list(phosphorylations)


def check_id(uniprot_id, seq, swiss_fasta):
    result = False
    path = os.path.join(swiss_fasta, uniprot_id + '.fasta')
    uni_seq = ''
    if os.path.exists(path):
        with open(path) as a:
            uni_seq = ''.join(a.read().splitlines()[1:])
    else:
        req = urllib2.Request("http://www.uniprot.org/uniprot/"
                              + uniprot_id + ".fasta")
        try:
            uni_seq = urllib2.urlopen(req).read()
        except urllib2.HTTPError and httplib.BadStatusLine:
            print "No entry with ID: {}".format(uniprot_id)
        else:
            uni_seq = ''.join(uni_seq.splitlines()[1:])
    if uni_seq == seq:
        result = True
    return result


def get_uniprot_txt(uniprot_id, swiss_txt):
    features = []
    go_terms = []
    path = os.path.join(swiss_txt,
                        uniprot_id + '.txt')
    uniprot_dat = []
    if os.path.exists(path):
        with open(path) as a:
            uniprot_dat = a.read().splitlines()
    else:
        try:
            req = urllib2.Request("http://www.uniprot.org/uniprot/{}.txt".format(
                uniprot_id))
            uniprot_dat = urllib2.urlopen(req).read().splitlines()
        except urllib2.HTTPError and httplib.BadStatusLine:
            print "Couldn't find a uniprot entry for id {}".format(uniprot_id)
    for lineI in uniprot_dat:
        if lineI.startswith('FT'):
            features += [lineI]
        elif lineI.startswith('DR   GO;'):
            go_terms += [lineI.split(';')[1].split(':')[1]]
    return {"features": features, "GO": go_terms}


def get_annotation_level(uni_features):
    levels_dict = {'269': 0, '314': 0, '353': 0, '315': 0, '316': 0, '270': 0,
                   '250': 1, '266': 1, '247': 1, '255': 1, '317': 1, '318': 1,
                   '319': 1, '320': 1, '321': 1, '245': 1,
                   '304': 2, '303': 2, '305': 2, '307': 2,
                   '501': 3}
    i = 0
    n = 3
    reading = True
    while reading and i < len(uni_features):
        if uni_features[i].startswith("FT          ") or i == 0:
            if "ECO:0000" in uni_features[i]:
                start = uni_features[i].index("ECO:0000") + 8
                eco_code = uni_features[i][start:start+3]
                if eco_code in levels_dict.keys():
                    n = levels_dict[eco_code]
                break
            i += 1
        else:
            reading = False
    return n


def get_ptm_type(line):
    ptm_found = True
    ptm = ""
    if "Phospho" in line and "MOD_RES" in line:
        ptm = 'phosph'
    elif "amide" in line and "MOD_RES":
        ptm = 'amids'
    elif "CARBOHYD" in line and "O-linked" in line:
        ptm = 'Oglycs'
    elif "CARBOHYD" in line and "N-linked" in line:
        ptm = 'Nglycs'
    elif "MOD_RES" in line and "acetyl" in line:
        ptm = 'acetyl'
    elif "MOD_RES" in line and "hydroxy" in line:
        ptm = 'hydrox'
    elif "MOD_RES" in line and "methyl" in line:
        ptm = 'meth'
    else:
        ptm_found = False
    return ptm, ptm_found


# UNIPROT
def find_ptm_sites(features):
    ptms = ["phosph", "Nglycs", "Oglycs", "amids", "hydrox",
            "meth", "acetyl"]
    ptms_dict = {i: [[] for j in range(4)] for i in ptms}
    for i, lineI in enumerate(features):
        if len(lineI.split()) > 3 and lineI.split()[3].isdigit():
            ptm, ptm_found = get_ptm_type(lineI)
            if ptm_found:
                if len(features) - i < 10:
                    end = len(features)
                else:
                    end = i+10
                n = get_annotation_level(features[i:end])
                position = int(lineI.split()[3])
                ptms_dict[ptm][n].append(position)
    return [ptms_dict['Oglycs'], ptms_dict['meth'], ptms_dict['hydrox'],
            ptms_dict['amids'], ptms_dict['Nglycs'], ptms_dict['acetyl'],
            ptms_dict['phosph']]


# filterOutOverlapping -> removes overlapping slims
# (the ones with lower probabilities are removed)
def filter_out_overlapping(lims, ids, probs):
    probs_with_indexes = [[i, probs[i]] for i in range(len(probs))]
    probs_with_indexes.sort(key=itemgetter(1))
    probs_with_indexes.reverse()
    new_lims = []
    new_ids = []
    new_probs = []
    for i in probs_with_indexes:
        goed = True
        start = lims[i[0]][0]
        end = lims[i[0]][1]
        for j in new_lims:
            if ((start >= j[0] and start <= j[1])
                    or (end >= j[0] and end <= j[1])):
                goed = False
        if goed:
            new_lims.append(lims[i[0]])
            new_ids.append(ids[i[0]])
            new_probs.append(probs[i[0]])
    return new_lims, new_ids, new_probs


def get_annotated_motifs(uniprotID):
    limits = []
    elms_ids = []
    probabilities = []
    # get annotated motifs first
    try:
        req = urllib2.Request("http://elm.eu.org/elms/browse_instances.gff"
                              "?q="+uniprotID)
        response = urllib2.urlopen(req)
        features = response.read().splitlines()

        for i in features:
            if 'sequence_feature' in i:
                start = int(i.split()[3])
                end = int(i.split()[4])
                elm_id = i.split()[8].split('=')[1]
                limits.append([start, end])
                probabilities.append(1)
                elms_ids.append(elm_id)
    except urllib2.HTTPError and httplib.BadStatusLine:
        print "can't get annotated motifs for {}".format(uniprotID)
    return [limits, elms_ids, probabilities]


def get_predicted_motifs(sequence, slims_all_classes, seq_go_terms):
    limits = []
    elms_ids = []
    probabilities = []
    try:
        data = urllib.urlencode({'sequence': sequence})
        url = "http://elm.eu.org/start_search/"
        req = urllib2.Request(url, data)
        response = urllib2.urlopen(req)
        features = response.read()
        features = features.splitlines()
        for line in features[1:]:
            entry = line.split()
            prob = 1
            if entry:
                slim_id = entry[0]
                try:
                    slim_go_terms = slims_all_classes[slim_id]["GO"]
                    if set(seq_go_terms).intersection(set(slim_go_terms)):
                        if entry[3] == "False":
                            prob = 1 + 1/math.log(
                                slims_all_classes[slim_id]["prob"], 10)
                            if prob > 0:
                                limits.append([int(entry[1]), int(entry[2])])
                                elms_ids.append(entry[0])
                                probabilities.append(prob)
                except KeyError:
                    print "Didn't find motif: {}".format(slim_id)
    except urllib2.HTTPError and httplib.BadStatusLine:
        print "can't get predicted motifs"
    return [limits, elms_ids, probabilities]


def search_elm(uniprotID, sequence, slims_all_classes, seq_go_terms):
    annotated = get_annotated_motifs(uniprotID)
    predicted = get_predicted_motifs(sequence, slims_all_classes, seq_go_terms)
    limits = annotated[0] + predicted[0]
    elms_ids = annotated[1] + predicted[1]
    probabilities = annotated[2] + predicted[2]
    limits, elms_ids, probabilities = filter_out_overlapping(limits,
                                                             elms_ids,
                                                             probabilities)
    return [limits, elms_ids, probabilities]


# adds new domains/motifs to the dictionary
# value is the encoded index (ab, ac, ...), key is the element's id
def add_elements_to_dict(newList, oldDict, dictkind):
    newDict = oldDict
    numbering = create_numbering()
    for i, elI in enumerate(newList):
        if len(elI) > 0:
            elID = elI.split('.')[0] if dictkind == "domains" else elI
            if elID not in newDict.keys():
                index = len(newDict.keys())
                encoded_index = numbering[index]
                newDict[elID] = encoded_index
    return newDict


# creates a list of codes for domains in myList
def get_codes(myDict, myList, mode):
    if mode == 'domains':
        result = [myDict[i.split('.')[0]] for i in myList]
    else:
        result = [myDict[i] for i in myList]
    return result


def encode_domains(seq, domains, domain_codes):
    for i, domI in enumerate(domains):
        start = domI[0]
        end = domI[1]
        domainsCode = domain_codes[i]
        j = 0
        k = 0
        while j < end+1 and k < len(seq):
            if k % 7 == 2 and j >= start and seq[k - 2] != '-':
                seq[k] = domainsCode[0]
                seq[k+1] = domainsCode[1]
            elif k % 7 == 0 and seq[k] != "-":
                j += 1
            k += 1
    return seq


def encode_predicted_phosph(seq, pred_phosph):
    # predictions are 1 based
    if pred_phosph:
        k = 0
        l = 0
        end = max(pred_phosph)
        while k < end + 1 and l < len(seq):
            if l % 7 == 4 and k in pred_phosph:
                seq[l] = "d"
            if l % 7 == 0 and seq[l] != "-":
                k += 1
            l += 1
    return seq


def encode_ptms(seq, ptms):
    code_table = [["Z", "a", "b", "c"],  # code table for PTMs
                  ["V", "W", "X", "Y"],
                  ["R", "S", "T", "U"],
                  ["J", "K", "L", "M"],
                  ["F", "G", "H", "I"],
                  ["B", "C", "D", "E"],
                  ["N", "O", "P", "Q"]]
    for i, ptmI in enumerate(ptms):
        for j, ptmModeJ in enumerate(ptmI):
            if ptmModeJ:
                k = 0
                l = 0
                end = max(ptmModeJ)
                while k < end + 1 and l < len(seq):
                    if l % 7 == 4 and k in ptmModeJ:
                        seq[l] = code_table[i][j]
                    if l % 7 == 0 and seq[l] != "-":
                        k += 1
                    l += 1
    return seq


def encode_slims(seq, slims, slim_codes):
    for i, slimI in enumerate(slims):
        if len(slimI) > 0:
            start = slimI[0]
            end = slimI[1]
            slimsCode = slim_codes[i]
            k = 0
            j = 0
            while j < end + 1 and k < len(seq):
                if k % 7 == 5 and j >= start and seq[k - 5] != '-':
                    seq[k] = slimsCode[0]
                    seq[k+1] = slimsCode[1]
                elif k % 7 == 0 and seq[k] != "-":
                    j += 1
                k += 1
    return seq


def sevenCharactersCode(results, myseq, domain_codes, slim_codes):
    """
    code: "012345" 0 - aa; 1 - nothing yet; 2 - domain, 3 - phosph; 4,5 - motif
    argument 'results' -> list([domains,
                                phosphorylations,
                                methods,
                                low complexity regions])
    """
    newseq = ""
    # first encode sequence with no features
    for i in myseq:
        newseq += i + "AAAAAA"
    # map results onto the sequence
    newseq = list(newseq)
    newseq = encode_domains(newseq, results[0], domain_codes)
    newseq = encode_predicted_phosph(newseq, results[4])
    newseq = encode_ptms(newseq, results[3])
    newseq = encode_slims(newseq, results[1], slim_codes)
    return ''.join(newseq)


# creates a tmp fastafile and returns a path to it
def tmp_fasta(seq_id, seq):  # pragma: no cover
    fasta_seq = ">{}\n{}\n".format(seq_id, seq)
    tmp_file = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
    with tmp_file as f:
        f.write(fasta_seq)
    return tmp_file.name


def elm_db(ELM_DB):
    with open(ELM_DB) as a:
        slims_all_classes_pre = a.read()
    return process_slims_all_classes(slims_all_classes_pre.splitlines())


def get_uniprot_data(fastafile, swiss_txt):
    result = dict({"seq_data": dict({}), "GO": set([])})
    for i in fastafile:
        if i.startswith('>'):
            seq_id = get_id(i.rstrip('\n'))
            seq_data = get_uniprot_txt(seq_id, swiss_txt)
            result["seq_data"][seq_id] = seq_data
            result["GO"] = result["GO"].union(seq_data["GO"])
    return result


# MAIN
def convert_to_7chars(filename, outname, ELM_DB, swiss_fasta, swiss_txt):
    slims_all_classes = elm_db(ELM_DB)
    seqFASTA = read_fasta(filename)
    domainsDictionary = dict()
    motifsDictionary = dict()
    motifProbsDict = dict()
    newfile = ""
    uniprot_data = get_uniprot_data(seqFASTA, swiss_txt)
    for i, seqI in enumerate(seqFASTA):
        if '>' not in seqI:
            seqI = seqI.rstrip("\n")
            seq_id = get_id(seqFASTA[i-1]).rstrip('\n')
            header = seqFASTA[i-1].rstrip('\n')
            ungapped = re.sub('-', '', seqI)
            print 'Annotating sequence: {}'.format(seq_id)
            tmp_filename = tmp_fasta(seq_id, seqI)
            [pfam, domains] = run_pfam_scan(tmp_filename)
            predicted_phosph = run_netphos(tmp_filename)
            if os.path.exists(tmp_filename):        # pragma: no cover
                os.remove(tmp_filename)
            if check_id(seq_id, ungapped, swiss_fasta):
                # uniprot_txt = get_uniprot_txt(seq_id)
                uniprot_results = find_ptm_sites(uniprot_data["seq_data"][seq_id]["features"])
                # elms - slims coordinates, motifs_ids, probs - probabilities
                [elm, motifs_ids, probs] = search_elm(seq_id, ungapped,
                                                      slims_all_classes,
                                                      uniprot_data["GO"])
                motifsDictionary = add_elements_to_dict(motifs_ids,
                                                        motifsDictionary,
                                                        'slims')
                motifs_codes = get_codes(motifsDictionary, motifs_ids, 'motifs')
                for i, indI in enumerate(motifs_codes):
                    motifProbsDict[indI] = probs[i]
            else:
                uniprot_results = []
                elm = []
                motifs_codes = []

            lc_regions = []
            domainsDictionary = add_elements_to_dict(domains,
                                                     domainsDictionary,
                                                     'domains')
            domains_codes = get_codes(domainsDictionary, domains, 'domains')
            resultsList = [pfam, elm, lc_regions,
                           uniprot_results, predicted_phosph]
            newfile += '{}\n{}\n'.format(header,
                                         sevenCharactersCode(resultsList,
                                                             seqI,
                                                             domains_codes,
                                                             motifs_codes))
    # write encoded fasta
    newfile += "## PROBABILITIES\n"
    for i in motifProbsDict:
        newfile += str(i)+' '+str(motifProbsDict[i])+'\n'

    out = open(outname, 'w')
    out.write(newfile)
    out.close()
    # write mapping of feature codes to their real names
    newfile = ''
    for i in motifsDictionary:
        newfile += 'motif_{} {} {}\n'.format(motifsDictionary[i], i,
                                             slims_all_classes[i]["regex"])
    for i in domainsDictionary:
        newfile += 'domain_{} {}\n'.format(domainsDictionary[i], i)
    out = open(filename.split('.')[0]+'.map', 'w')
    out.write(newfile)
    out.close()


if __name__ == "__main__":
    logging.basicConfig(level=logging.WARNING)
    log = logging.getLogger('convert')
    script_dir = os.path.dirname(os.path.realpath(__file__))
    src_dir = '/'.join(script_dir.split('/')[:-1])
    ELM_DB = src_dir + '/data/' + 'elm_complete.txt'
    parser = argparse.ArgumentParser(description='Convert fasta to a KMAD'
                                                 + ' compatible format')
    parser.add_argument('input_filename')
    parser.add_argument('output_filename')
    parser.add_argument('--elm_db', nargs='?', default=ELM_DB)
    parser.add_argument('--swiss_fasta', nargs='?', default="dummy")
    parser.add_argument('--swiss_txt', nargs='?', default="dummy")
    args = parser.parse_args()
    if not os.path.exists(args.elm_db):
        try:
            raise Exception("ELM directory not found (path: {})".format(
                args.elm_db))
        except Exception, err:
            log.exception(err)
            sys.exit(3)
    convert_to_7chars(args.input_filename, args.output_filename, args.elm_db,
                      args.swiss_fasta, args.swiss_txt)
