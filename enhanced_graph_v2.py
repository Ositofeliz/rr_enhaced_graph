#!/usr/bin/env python3

###########################################################################################

# Programme d'améliorations des graphes de contigs par l'exploitation des régions répétées

# Auteurs : Quentin Delorme, Rémy Costa (qdelorme@yahoo.fr, remy.costa@live.fr)

# Réalisé sous l'encadrement de Annie Chateau, LIRMM, équipe MAB

###########################################################################################

import os, sys, re, argparse, pprint, random
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from Bio import SeqIO


# ARGUMENTS

parser = argparse.ArgumentParser(description="Clustering contigs by Repeated Regions (RR)")

# Blast file output
parser.add_argument("-b", "--blast", dest = "blast_out", type = str, help = "Repeated regions blast output file")

# Contigs fasta file
parser.add_argument("-seq", "--seq_fa", dest = "contig_file", type = str, help = "Contigs FASTA file")

# Repeated regions sequences reference file
parser.add_argument("-ref", "--refbase", dest = "reference", type = str, help = "Path to reference database")

# Paired-end graph dot file
parser.add_argument("-dot", "--dotfile", dest = "dotfile", type = str, help = "dot file")

args = parser.parse_args()

# Catching arguments
blast_out = args.blast_out
contig_file = args.contig_file
ref = args.reference
dot_file = args.dotfile

pp = pprint.PrettyPrinter(indent=4)


###########################################################################################
#                       FLAGGING CONTIGS WITH NON-FACTORISED RR                           #
###########################################################################################

# Creating transposable elements dictionnary 
rr_info = {}    # key : contig_id, value : RR type + positions in the contig

with open(blast_out) as blast_frame :
    for line in blast_frame :
        id = line.split()[0]
        rr_info[id] = line.split()[1] + '\t' + line.split()[6] + '-' + line.split()[7]

# Reading contigs file
with open(contig_file) as contig_frame :
    # Edition fasta file by adding contig informations in headers of the fasta file
    target_file = open("./RR_annotated_seq.fa", 'w')    # Annotated output fasta file
    print("Editing contig informations")

    contig_fasta = open(contig_file, 'r')
    for rec in SeqIO.parse(contig_fasta, "fasta") :
        if rec.id in rr_info :
            target_file.write('>' + rec.id + '\t' + rr_info[rec.id] + '\t' + str(len(rec)) + '\t' + '\n' + str(rec.seq) + '\n')
        else :
            target_file.write('>' + rec.id + '\t' + str(len(rec)) + '\n' + str(rec.seq) + '\n')

# Afficher rr_info
#pp.pprint(rr_info)
#sys.exit(0)
###########################################################################################
#                DATABASE CONCATENATION AND CATCHING REPEATS FAMILY NAMES                 #
###########################################################################################

count_entry = 0
fam_dic = {}

with open(ref) as ref_frame :
    # Counting number of entries in RR families
    for line in ref_frame :
        if line[0] == '>' :
            count_entry += 1

# Editing IDs
    rand_id = random.sample(range(count_entry * 2), k = count_entry)
    out_target = open("ref_repreg.fasta", 'w')
    ite_entry = 0
    # Creating repeated regions families dictionnary
    ref_frame.seek(0) # replace le curseur au début du fichier
    for line in ref_frame :
        if line[0] == '>' and len(line.split('\t')) > 1 :
            rep_id = line.split('\t')[0].replace('>', '')
            fam_name = line.split('\t')[1]
            fam_dic[rep_id] = fam_name + '\t' + str(rand_id[ite_entry])
            out_target.write(line.replace('\n', '') + '\t' + str(rand_id[ite_entry]) + '\n')
            ite_entry += 1
        else :
            out_target.write(line)


#pp.pprint(fam_dic)
#sys.exit(0)

###########################################################################################
#                                  REPEATS FACORIZATION                                   #  
###########################################################################################

print("Factorizing repeats")

# Counting repeats
rep_dic = {}

with open("./RR_annotated_seq.fa") as annorr_frame:
    # Counting repeats
    for line in annorr_frame :
        if line[0] == '>' and len(line.split('\t')) > 2 :
            repeat = line.split('\t')[1].replace('\n', '')
            if repeat not in rep_dic :
                rep_dic[repeat] = 1
            else :
                rep_dic[repeat] += 1
        else :
            continue

    # Sequence file by repeats with factorized families
    target_file = open("./factorized_contig.fa", 'w')
    fact_rep = {}

    annorr_frame.seek(0)
    for line in annorr_frame :
        if line[0] == '>' and len(line.split('\t')) > 2 :
            rep = line.split('\t')[1].replace('\n', '')
            not_seen = True
            if rep in fam_dic :
                target_file.write(line.split('\t')[0] + '\t' + fam_dic[rep] + '\t' + line.split('\t')[2] + '\t' + str(line.split('\t')[3]) + '\n')
                not_seen = False
                if not_seen :
                    target_file.write(line)
        else : 
            target_file.write(line)

#pp.pprint(rep_dic)
#sys.exit(0)

###########################################################################################
#                                 ANALYZING FACOTRIZED FILE                               #  
###########################################################################################

from repeats import Repeat

# Creating a list containing Repeat objects
repeat_obj = []

with open ("./factorized_contig.fa", 'r') as fact_frame :
    # Analyzing file
    for line in fact_frame :
        if line[0] == '>' :
            if len(line.split('\t')) > 2 :
                if line.split('\t')[-1] != '\n' :
                    contig_name = line.split('\t')[0].replace('>', '')
                    rep_name = line.split('\t')[1]
                    rep_id = line.split('\t')[2]
                    start = line.split('\t')[3].replace('\n', '').split('-')[0]
                    end = line.split('\t')[3].replace('\n', '').split('-')[1]
                    contig_len = str(line.split('\t')[4]).replace('\n', '')
                    seq = "none"
                    if '/' in rep_name :
                        valid_rname = rep_name.replace('/', '_')
                        repeat_obj.append(Repeat(valid_rname, rep_id, contig_name, contig_len, start, end, seq))
                    else :
                        repeat_obj.append(Repeat(rep_name, rep_id, contig_name, contig_len, start, end, seq))
            else :
                contig_name = line.replace('>', '').replace('\n', '').split('\t')[0]
                rep_name = "none"
                rep_id = "none"
                start = "none"
                end = "none"
                seq = "none"
                contig_len = str(line.split('\t')[1].replace('\n', ''))
                repeat_obj.append(Repeat(rep_name, rep_id, contig_name, contig_len, start, end, seq))
        else :
            repeat_obj[-1].set_contig_seq(line.replace('\n', ''))

cluster_rep = {}
for elt in repeat_obj :
    if elt.get_rep_name() not in cluster_rep :
        cluster_rep[elt.get_rep_name()] = [elt]
    else :
        cluster_rep[elt.get_rep_name()].append(elt)

#pp.pprint(repeat_obj)
#pp.pprint(cluster_rep)
#sys.exit(0)

###########################################################################################
#                                    GENERATING RR GRAPHS                                 #
###########################################################################################

# Creating list containing RR edges
#rr_graph_edges = []
count_rr_graph_edges = 0
ite_rr_edges = 0

for famrep in cluster_rep.keys() :   # For each family cluster
    if famrep == "none" :
        continue
    if not os.path.isdir('./rr_graph') :
        print("Generating files containing repeat cluster graphs")
        os.makedirs("./rr_graph")
    print("Generating cluster repeat family graph " + famrep)
    graph_file = open("./rr_graph/{}_graph.dot".format(famrep), 'w')
    graph_file.write("digraph {}_graph".format(famrep) + '{' + '\n')
    for i in range(0, (len(cluster_rep[famrep])-1),1) : # Contig n
        for j in range(i+1, len(cluster_rep[famrep]), 1) : # Contig n+1
            #print(cluster_rep[famrep][j].rep_id)
	    #print(cluster_rep[famrep][i].rep_id)
            if cluster_rep[famrep][j].rep_id == cluster_rep[famrep][i].rep_id :
                if cluster_rep[famrep][j].pos == "ext" and cluster_rep[famrep][i].pos == "ext" :
                    if int(cluster_rep[famrep][i].start) == 1 or int(cluster_rep[famrep][i].end) == 1 :
                        if int(cluster_rep[famrep][j].start) == 1 or int(cluster_rep[famrep][j].end) == 1 :
                            edge = cluster_rep[famrep][i].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][j].contig_name.replace(" ", '')
                            count_rr_graph_edges += 1
                            graph_file.write('"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "- to +" + '"' + "]" + ';' + '\n')
                        if int(cluster_rep[famrep][j].start) == int(cluster_rep[famrep][j].contig_len) or int(cluster_rep[famrep][j].end) == int(cluster_rep[famrep][j].contig_len) :
                            edge = cluster_rep[famrep][j].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][i].contig_name.replace(" ", '')
                            count_rr_graph_edges += 1
                            graph_file.write('"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "+ to +" + '"' + "]" + ';' + '\n')
                    if int(cluster_rep[famrep][j].start) == 1 or int(cluster_rep[famrep][j].end) == 1 :
                        if int(cluster_rep[famrep][i].start) == 1 or int(cluster_rep[famrep][i].end) == 1 :
                            edge = cluster_rep[famrep][j].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][i].contig_name.replace(" ", '')
                            count_rr_graph_edges += 1
                            graph_file.write('"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "- to +" + '"' + "]" + ';' + '\n')
                        if int(cluster_rep[famrep][i].start) == int(cluster_rep[famrep][i].contig_len) or int(cluster_rep[famrep][i].end) == int(cluster_rep[famrep][i].contig_len) :
                            edge = cluster_rep[famrep][i].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][j].contig_name.replace(" ", '')
                            count_rr_graph_edges += 1
                            graph_file.write('"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "+ to +" + '"' + "]" + ';' + '\n')
                    if int(cluster_rep[famrep][i].start) == int(cluster_rep[famrep][i].contig_len) or int(cluster_rep[famrep][i].end) == int(cluster_rep[famrep][i].contig_len) :
                        if int(cluster_rep[famrep][j].start) == 1 or int(cluster_rep[famrep][j].end) == 1 :
                            edge = cluster_rep[famrep][i].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][j].contig_name.replace(" ", '')
                            count_rr_graph_edges += 1
                            graph_file.write('"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "+ to +" + '"' + "]" + ';' + '\n')
                        if int(cluster_rep[famrep][j].start) == int(cluster_rep[famrep][j].contig_len) or int(cluster_rep[famrep][j].end) == int(cluster_rep[famrep][j].contig_len) :
                            edge = cluster_rep[famrep][i].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][j].contig_name.replace(" ", '')
                            count_rr_graph_edges += 1
                            graph_file.write('"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "+ to -" + '"' + "]" + ';' + '\n')
                    if int(cluster_rep[famrep][j].start) == int(cluster_rep[famrep][j].contig_len) or int(cluster_rep[famrep][j].end) == int(cluster_rep[famrep][j].contig_len) :
                        if int(cluster_rep[famrep][i].start) == 1 or int(cluster_rep[famrep][i].end) == 1 :
                            edge = cluster_rep[famrep][j].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][i].contig_name.replace(" ", '')
                            count_rr_graph_edges += 1
                            graph_file.write('"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "+ to +" + '"' + "]" + ';' + '\n')
                        if int(cluster_rep[famrep][i].start) == int(cluster_rep[famrep][i].contig_len) or int(cluster_rep[famrep][i].end) == int(cluster_rep[famrep][i].contig_len) :
                            edge = cluster_rep[famrep][j].contig_name.replace(" ", '') + '\t' + cluster_rep[famrep][i].contig_name.replace(" ", '')                            
                            count_rr_graph_edges +=1
                            graph_file.write('"' + cluster_rep[famrep][j].contig_name.replace(" ", '') + '"' + '->' + '"' + cluster_rep[famrep][i].contig_name.replace(" ", '') + '"' + "[cname=" + '"' + str(cluster_rep[famrep][j].rep_id) + '"' + ",label=" + '"' + "+ to -" + '"' + "]" + ';' + '\n')

    graph_file.write("} \n")
    
#pp.pprint(rr_graph_edges)
#sys.exit()

###########################################################################################
#                                    GRAPHS CROSSING                                      #
###########################################################################################

# Creating dictionnary containing contigs' edges and vertices
dot_contigs_desc = {}
# Searching for descriptive lines of contigs' edges 
contig_edges_search = re.compile("([0-9]*)--([0-9]*).*\"(.*)\".*")

# Creating dictionnary containing intercontigs edges (normal and reverse)
dot_inter_desc = {}
inter_edge_search = re.compile("([0-9]*) -- ([0-9]*).*label=([0-9]*).*")

with open(dot_file) as dot_frame :
    for line in dot_frame :
        pattern_search = contig_edges_search.search(line)
        inter_search = inter_edge_search.search(line)
        if pattern_search :
            vertex_start = str(pattern_search.group(1))
            vertex_end = str(pattern_search.group(2))
            contig_id = pattern_search.group(3)
            dot_contigs_desc[vertex_start] = contig_id
            dot_contigs_desc[vertex_end] = contig_id
        if inter_search :
            emitter = str(inter_search.group(1))
            receiver = str(inter_search.group(2))
            support = str(inter_search.group(3))
            edge = dot_contigs_desc[emitter] + '\t' + dot_contigs_desc[receiver]
            edge_alter = dot_contigs_desc[receiver] + '\t' + dot_contigs_desc[emitter]
            dot_inter_desc[edge] = support

    ite_common_edges = 0 
    common_edge = {}
    uncommon_edge = {}

    print("Number of edges in Paired-End Graph : " + str(len(dot_inter_desc)))
    print("Number of edges in Repeat Graphs : " + str(count_rr_graph_edges))

    # PARTIE MODIFIEE POUR UTILISATION MEMOIRE
    rr_graphs = os.listdir("rr_graph/")
    for rr_file in rr_graphs :
        with open("rr_graph/{}".format(rr_file)) as current_graph :
            for line in current_graph :
                res = re.search("(NODE.*)\"->\"(NODE.*)\"\[", line)
                if res : 
                    rr_edge = res.group(1) + '\t' + res.group(2)
                    if rr_edge in dot_inter_desc :
                        ite_common_edges += 1
                        common_edge[rr_edge] = dot_inter_desc[edge]
                    else :
                        uncommon_edge[rr_edge] = "0"

    #pp.print(common_dege)
    print("Number of shared edges (dot vs RR) : " + str(len(common_edge)))

    stat_file = open("RRPE_stats.csv" , 'w')
    stat_file.write("afferent" + '\t' + "efferent" + '\t' +  " PE score" + '\t' + "repeat type" + '+\t' + "size of cluster" + '\n')
    common_with_clustlen = {}

    # AMELIORATION DU POIDS DES ARETES - GITAN ++
    for edge_co in common_edge :
        score_support = str(common_edge[edge_co])
        emitter_contig = edge_co.split('\t')[0]
        receiver_contig = edge_co.split('\t')[1]
        for rep_clust in cluster_rep :
            for contig_obj in cluster_rep[rep_clust] :
                if emitter_contig == contig_obj.contig_name :
                    stat_file.write(emitter_contig + '\t' + receiver_contig + '\t' + score_support + '\t' + rep_clust + '\t' + str(len(cluster_rep[rep_clust])) + '\n')
                    eucli_score = len(cluster_rep[rep_clust]) // 100
                    common_with_clustlen[edge_co] = eucli_score
                    
    # Edges linking two contigs carrying RR
    inter_to_RR = {}
    for elt in dot_inter_desc :
        emitter = elt.split('\t')[0]
        receiver = elt.split('\t')[1]
        if emitter in rr_info and receiver in rr_info :
            inter_to_RR[elt] = dot_inter_desc[elt]
    print("Two contigs carrying RR : " + str(len(inter_to_RR)))

    # Eliminating edges from dot file because of their different RRs :
    eliminated_edge = {}
    for elt in inter_to_RR :
        emitter = elt.split('\t')[0]
        receiver = elt.split('\t')[1]
        if rr_info[emitter].split('\t')[0] != rr_info[receiver].split('\t')[0] :
            eliminated_edge[elt] = inter_to_RR[elt]
    print("Two contigs carrying RRs eliminated : " + str(len(eliminated_edge)))

    # Eliminating edges from dot file linking a contig carrying RR and a non-carrying RR contig
    inter_RRtonone = {}
    for elt in dot_inter_desc :
        emitter = elt.split('\t')[0]
        receiver = elt.split('\t')[1]
        if emitter in rr_info and receiver not in rr_info :
            inter_RRtonone[elt] = dot_inter_desc[elt]
    print("One carrying-RR contigs : " + str(len(inter_RRtonone)))

    # Eliminating edges from dot file linking a non-carrying RR contig to a carrying RR contig
    elim_nonetoRR = {}
    for elt in repeat_obj :
        con_name = elt.contig_name
        for elem in inter_RRtonone :
            emitter = elem.split('\t')[0]
            receiver = elem.split('\t')[1]
            if con_name == emitter :
                if elt.pos == "ext" and int(elt.start) == int(elt.contig_len) :
                    elim_nonetoRR[elem] = inter_RRtonone[elem]
                if elt.pos == "ext" and int(elt.end) == int(elt.contig_len) :
                    elim_nonetoRR[elem] = inter_RRtonone[elem]
            if con_name == receiver :
                if elt.pos == "ext" and int(elt.start) == 1 :
                    elim_nonetoRR[elem] = inter_RRtonone[elem]
                if elt.pos == "ext" and int(elt.end) == 1 :
                    elim_nonetoRR[elem] = inter_RRtonone[elem]
    print("One carrying-RR contigs eliminated : " + str(len(elim_nonetoRR)))

    # Generating shared contigs fasta file
    repeat_obj_dic = {}
    for elt in repeat_obj :
        con_name = elt.contig_name
        repeat_obj_dic[con_name] = elt.contig_seq

    common_arete = open("Contigs_from_common_edges.txt" , 'w')
    print("Generating file containing info on shared contig edges : ./Contigs_from_common_edges.txt")
    for elem in common_edge :
        emitter = elem.split('\t')[0]
        receiver = elem.split('\t')[1]
        support = common_edge[elem]
        if emitter in repeat_obj_dic and receiver in repeat_obj_dic :
            common_arete.write(emitter + '\t' + receiver + '\t' + support + '\t' + '\n')

    RR2RR_edge = open("Contigs_from_RvR_edge.txt" , 'w')
    print("Generating file containing info on RvR contig edges : ./Contigs_from_RvR_edge.txt")
    for elem in eliminated_edge:
        emitter = elem.split('\t')[0]
        receiver = elem.split('\t')[1]
        support = eliminated_edge[elem]
        if emitter in repeat_obj_dic and receiver in repeat_obj_dic :
            RR2RR_edge.write(emitter + '\t' + receiver + '\t' + support + '\n')

    N2RR_edge = open("Contigs_from_NvR_edge.txt" , 'w')
    print("Generating file containing info on NvR contig edges : ./Contigs_from_NvR_edge.txt")
    for elem in elim_nonetoRR:
        emitter = elem.split('\t')[0]
        receiver = elem.split('\t')[1]
        support = elim_nonetoRR[elem]
        if emitter in repeat_obj_dic and receiver in repeat_obj_dic :
            N2RR_edge.write(emitter + '\t' + receiver + '\t' + support + '\n')

    print("Generating enhanced paired-end graph : ./enhanced_PE_graph.dot")
    enhanced_dot = open("enhanced_PE_graph.dot" , 'w')
    dot_frame.seek(0)
    for line in dot_frame :
        inter_search = inter_edge_search.search(line)
        if inter_search :
            emitter = str(inter_search.group(1))
            receiver = str(inter_search.group(2))
            support = str(inter_search.group(3))
            edge = dot_contigs_desc[emitter] + '\t' + dot_contigs_desc[receiver]
            anti_edge = dot_contigs_desc[receiver] + '\t' + dot_contigs_desc[emitter]
            if edge in common_with_clustlen or anti_edge in common_with_clustlen :
                if edge in common_with_clustlen :
                    new_score = common_with_clustlen[edge] + int(support)
                    p = re.compile("label=[0-9]*")
                    new_line = p.sub("label="+str(new_score), line).replace("];", ",color=gold];")
                    enhanced_dot.write(new_line)
                if anti_edge in common_with_clustlen:
                    new_score = common_with_clustlen[anti_edge] + int(support)
                    p = re.compile("label=[0-9]*")
                    new_line = p.sub("label="+str(new_score), line).replace("];", ",color=gold];")
                    enhanced_dot.write(new_line)
            else :
                if edge not in eliminated_edge and anti_edge not in eliminated_edge:
                    if anti_edge not in elim_nonetoRR and anti_edge not in elim_nonetoRR:
                        enhanced_dot.write(line)
        else:
            enhanced_dot.write(line)

#pp.pprint(dot_contigs_desc)
#pp.pprint(dot_inter_desc)
