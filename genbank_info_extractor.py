#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This algorithm extracts information from files in GenBank(.gb) format,
generating one file in tabular format(TSV) and another in FASTA format with
the desired parameters.'''

__author__ = "Ary Rivillas"
__copyright__ = "Copyright 2020, Biotecnologia Microbiana Research Group,\
     Universidad Nacional de Colombia - Sede Medellín"
__credits__ = ["Pablo Gutierrez", "Daniel Tejada", "Andrea Restrepo",\
     "Susana Giraldo"]
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Ary Rivillas"
__email__ = "amrivillast@unal.edu.co"
__status__ = "Development"
__date__ = '2020/06/5'

import re 
import os

'''The file to be processed must be in the script folder. Modify the
following variable and enter the corresponding name'''
##############################################################################
############################### MODIFY ME ####################################
##############################################################################

file_base = "viruses_landplants_refseqs.gb"

##############################################################################
##############################################################################
##############################################################################

dir_path = os.path.dirname(os.path.realpath(__file__))
path_to_file = os.path.join(dir_path, file_base)

with open (path_to_file) as f:           
    
    """ The lists necessary to save the information of interest are created. 
       Also counters for cycles are created """
    length = []
    code = []
    sequence_l = []       
    sequence_temp = ""
    description = []
    isolate_l = []
    host_l = []
    country_l = []
    segment_l = []
    gc = []
    organism_l = []
    keyword_l = []
    family_l = []
    genus_l = []
    pudmed_l =[]        
    z = 0
    w = 0
            
    """The file is scrolled line by line, listing each line
    """
    for i, line in enumerate(f): 
        
        """The LOCUS pattern is searched, if found, the current line is
           saved in locus_hit
        """
        locus = '^LOCUS'
        locus_hit = re.match(locus, line)
        
        """If the regular expression is found, the third column of the line
           is selected and saved in a list
        """
        if locus_hit:
                   
            length.append((line.split()[2]))

        version = '^VERSION'
        version_hit = re.match(version, line)
        
        if version_hit:

            code.append((line.split()[1]))

        keyword = '^KEYWORDS'
        keyword_hit = re.match(keyword, line)

        if keyword_hit:
            
            line_mod1_1 = line.replace("KEYWORDS    ", "", 1)
            line_mod1_2 = line_mod1_1.replace(".", "", 1)
            line_mod1_3 = line_mod1_2.strip()
            
            """It is identified if the sequence to be processed is reference
               or unverified.
            """
            refs = re.findall(r'[R]+[e]+[f]+[S]+[e]+[q]+', line_mod1_3)
            unv = re.findall(r'[U]+[N]+[V]+[E]+[R]+[I]+[F]+[I]+[E]+[D]+',\
                             line_mod1_3)
            
            if (len(refs)) >= 1:
                            
                keyword_l.append("R")
            
            if (len(unv)) >= 1:
                            
                keyword_l.append("U")
        
        organism = '^  ORGANISM'
        organism_hit = re.match(organism, line)

        if organism_hit:

            line_mod2_1 = line.replace("ORGANISM  ", "", 1)
            line_mod2_2 = line_mod2_1.strip()
            organism_l.append(line_mod2_2)
            
            """This section is used to find the family and gender
            """
            next_line1_1 = next(f)                        # Go to the next line
            next_line1_2 = next_line1_1.strip()       
            line_find_1 = next_line1_2.split(';')
            next_line1_3 = next(f)            
            next_line1_4 = next_line1_3.strip()
            line_fine_2 = next_line1_4.split(';')
            line_find_f = line_fine_2 + line_find_1         # Concatenate lists

            for i in line_find_f:
                
                if 'viridae' in i and (len(code) - len(family_l) == 1):
                    
                    i_mod = i.replace(".","").strip()
                    family_l.append(i)
                    
                if 'virus' in i and (len(code) - len(genus_l) == 1):
                    
                    i_mod = i.replace(".","").strip()
                    genus_l.append(i_mod)
                    
            if (len(code) - len(family_l) == 1):
                
                family_l.append('nd')
                
            if (len(code) - len(genus_l) == 1):
                
                genus_l.append('nd')
        
        pubmed = '^   PUBMED'
        pubmed_hit = re.match(pubmed, line)

        if pubmed_hit:
            
            line_mod3_1 = line.strip()
            line_mod3_2 = line_mod3_1.replace("PUBMED", "")
            
            if (len(code) - len(pudmed_l) == 1):
                
                pudmed_l.append(line_mod3_2.strip())

        if (len(code) - len(pudmed_l) == 2):
            
            pudmed_l.append('nd')

        features = '^FEATURES'
        features_hit = re.match(features, line)
            
        if features_hit:
            
            """This cycle allows you to advance line by line until you find the
               end, that is, GENE           
            """
            while w == 0:
                
                next_line2_1 = next(f)
                stop_features = '^ORIGIN'
                stop_features_hit = re.match(stop_features, next_line2_1)
     
                if stop_features_hit:
                    
                    sequence_temp = ""
                    
                    while z == 0:
                                               
                        next_line2_2 = next(f)                         
                        stop_sequence_l = '^//'
                        stop_sequence_l_hit = re.match(stop_sequence_l, \
                                                       next_line2_2)
                                                
                        if stop_sequence_l_hit:
                            
                            z = 1
         
                        else:
                            
                            line_mod4_1 = (next_line2_2.split()[0])
                            line_mod4_2 = next_line2_2.replace(line_mod4_1,"",1)
                            line_mod4_3 = line_mod4_2.strip()
                            line_mod4_4 = line_mod4_3.replace(" ", "")
                            sequence_temp = sequence_temp + line_mod4_4
              
                    sequence_l.append(sequence_temp)
                    sequence_temp = ""
                    z = 0      
                    w = 1
                    
                else:

                    segment = '^                     /segment="'
                    segment_hit = re.match(segment, next_line2_1)
                    
                    isol = '^                     /isolate="'
                    isol_s_hit = re.match(isol, next_line2_1)
                    
                    host = '^                     /host="'
                    host_hit = re.match(host, next_line2_1)
                    
                    country = '^                     /country="'
                    country_hit = re.match(country, next_line2_1)
                    
                    if segment_hit:                             
                          
                        segment = next_line2_1.replace('/segment="',"",1)
                        segment_mod = segment.replace('"',"")
                        segment_mod2 = segment_mod.strip()
                        segment_l.append(segment_mod2)              
  
                    if isol_s_hit:                             
                          
                        isolate = next_line2_1.replace('/isolate="',"",1)
                        isolate_mod = isolate.replace('"',"")
                        isolate_mod2 = isolate_mod.strip()
                        isolate_l.append(isolate_mod2)
                                                                       
                    if host_hit:                             
                          
                        host = next_line2_1.replace('/host="',"",1)
                        host_mod = host.replace('"',"")
                        host_mod2 = host_mod.strip()
                        host_l.append(host_mod2)
                                                      
                    if country_hit:                             
                          
                        country = next_line2_1.replace('/country="',"",1)
                        country_mod = country.replace('"',"")
                        country_mod2 = country_mod.strip()
                        country_l.append(country_mod2)                       
                       
            w = 0
            
            if (len(code) - len(organism_l) == 1):
                
                organism_l.append("nd")
            
            if (len(code) - len(segment_l) == 1):
                
                segment_l.append("nd")                
                        
            if (len(code) - len(isolate_l) == 1):
                
                isolate_l.append("nd")            
                        
            if (len(code) - len(host_l) == 1):
                
                host_l.append("nd")            
                        
            if (len(code) - len(country_l) == 1):
                
                country_l.append("nd")
  
"""With this cycle each element of the sequence list is traversed in order to
   calculate the percentage guanine-cytosine present """
for i in sequence_l:
        
    a_s = i.count("a")             
    t_s = i.count("t")
    g_s = i.count("g")
    c_s = i.count("c")         
    t_gc = g_s + c_s
    t_gcat = g_s + c_s + t_s + a_s         
    per_gc = ((t_gc)/(t_gcat))*100
    gc.append(round(per_gc, 2))

"""Since all the generated lists have the same size, the elements are
   correlated by position so it is possible to create several dictionaries
   that relate the code to a certain characteristic """
dict_code_code = dict(zip(code, code))
dict_code_host = dict(zip(code, host_l))
dict_code_length = dict(zip(code, length))
dict_code_isolate = dict(zip(code, isolate_l))
dict_code_country = dict(zip(code, country_l))
dict_code_organism = dict(zip(code, organism_l))
dict_code_sequence = dict(zip(code, sequence_l))
dict_code_gc = dict(zip(code, gc))
dict_code_segment = dict(zip(code, segment_l))
dict_code_status = dict(zip(code, keyword_l))
dict_code_family = dict(zip(code, family_l))
dict_code_genus = dict(zip(code, genus_l))
dict_code_pudmed = dict(zip(code, pudmed_l))

file_base_tsv = "database_viral.tsv"
path_to_file_tsv = os.path.join(dir_path, file_base_tsv)

""" A tabulated file is created with the characteristics of interest in a 
   specific order. """
with open(path_to_file_tsv, "a+") as tsv: 
    
    for i in code:

        length_gen = dict_code_length.get(i)
        organism_gen  = dict_code_organism.get(i)
        code_gen  = dict_code_code.get(i)
        host_gen  = dict_code_host.get(i)
        isolate_gen = dict_code_isolate.get(i)
        country_gen = dict_code_country.get(i)
        seq_gen = dict_code_sequence.get(i)
        gc_gen = dict_code_gc.get(i)
        segment_gen = dict_code_segment.get(i)
        status_gen = dict_code_status.get(i)        
        family_gen = dict_code_family.get(i)
        genus_gen = dict_code_genus.get(i)
        pudmed_gen = dict_code_pudmed.get(i)
                
        print(str(code_gen) + '\t' + str(organism_gen) + '\t'\
              + str(family_gen) + '\t' + str(genus_gen) + '\t'\
              + str(pudmed_gen) + '\t' + str(length_gen) + '\t'\
              + str(host_gen) + '\t' + str(isolate_gen) + '\t'\
              + str(country_gen) + '\t' + str(segment_gen) + '\t'\
              + str(status_gen) + '\t' + str(gc_gen) + '\t'\
              + str(seq_gen), file = tsv)    

file_base_fasta = "database_viral.fasta"
path_to_file_fasta = os.path.join(dir_path, file_base_fasta)

""" A file is created in fasta format with the characteristics of interest in
   a specific order """
with open(path_to_file_fasta,"a+") as fasta: 
    
    for i in code:

        length_gen = dict_code_length.get(i)
        organism_gen  = dict_code_organism.get(i)
        code_gen  = dict_code_code.get(i)
        host_gen  = dict_code_host.get(i)
        isolate_gen = dict_code_isolate.get(i)
        country_gen = dict_code_country.get(i)
        seq_gen = dict_code_sequence.get(i)
        gc_gen = dict_code_gc.get(i)
        segment_gen = dict_code_segment.get(i)
        status_gen = dict_code_status.get(i)        
        family_gen = dict_code_family.get(i)
        genus_gen = dict_code_genus.get(i)
        pudmed_gen = dict_code_pudmed.get(i)
                
        print(">" + str(code_gen) + '\t' + str(organism_gen) + '\t'\
              + str(family_gen) + '\t' + str(genus_gen) + '\t'\
              + str(pudmed_gen) + '\t' + str(length_gen) + '\t'\
              + str(host_gen) + '\t' + str(isolate_gen) + '\t'\
              + str(country_gen) + '\t' + str(segment_gen) + '\t'\
              + str(status_gen) + '\t' + str(gc_gen), file = fasta)
        print(str(seq_gen), file = fasta)
                