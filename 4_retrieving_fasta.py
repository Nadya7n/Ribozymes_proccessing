#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 06.06.2021
#@author: Nadezhda Sorokina
#@contact: nadya7n@bk.ru

import os
from tqdm import tqdm
from Bio import SeqIO
from typing import Dict


def retrieve_ribo_coordinates(filtered_input_file: str, output_coord: str) -> None:
    with open(filtered_input_file) as fh:
        with open(output_coord, "a+") as fw:
            
            for record in SeqIO.parse(fh, "fasta"):
                header, info = record.description.split("\t")[0].split(":"), record.seq.split()
                id_org = header[0].split(">")[1]
                id_family = header[2]
                start, end, strand = info[9], info[10], info[11]

                
                fw.write(f"{id_org}:{start}:{end}:{strand}:{id_family}\n")


def parse_coord_from_ribo(filtered_coord_file: str, dct_with_coord_ribo=dict()) -> Dict[str, str]:
    with open(filtered_coord_file) as fh:
        for line in fh:
            line = line.strip().split(":")
            id_seq_from_ribo = line[0]
            coord_with_strand = line[1:5]
            if id_seq_from_ribo not in dct_with_coord_ribo:
                dct_with_coord_ribo[id_seq_from_ribo] = [coord_with_strand,]
            else:
                dct_with_coord_ribo[id_seq_from_ribo].append(coord_with_strand)
    return dct_with_coord_ribo


def write_ribo_fasta(data: str, output_file: str) -> None:
    with open(output_file, "a+") as fw:
        seq_id, start, end, strand ,ribo_id, seq_ribo = data
        fw.write(f">{seq_id}:{start}:{end}:{strand}:{ribo_id}\n")
        fw.write(f"{seq_ribo}\n")
        
        
def parse_fasta(filtered_coord_file: str, path_to_fasta_dir: str, path_to_output: str) -> None:
    dct_with_coord_ribo = parse_coord_from_ribo(filtered_coord_file)

    for id_g, fasta_file in enumerate(tqdm(os.listdir(path_to_fasta_dir))):
        fasta_file = os.path.join(path_to_fasta_dir, fasta_file)

        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            seq_id, seq = seq_record.id, seq_record.seq
            reverse_seq = seq.reverse_complement()
            
            if seq_id in dct_with_coord_ribo:
                for item in dct_with_coord_ribo[seq_id]:
                    start, end, strand, ribo_id = item
                    start, end = min(start,end), max(start, end)
                    seq_ribo = seq[int(start):int(end)] if strand == "+" else reverse_seq[int(start):int(end)]
                    write_ribo_fasta((seq_id, start, end, strand ,ribo_id, seq_ribo), path_to_output)

                        
if __name__ == "__main__":                       
    path_to_fasta_dir = "/path_to_fasta_dir/"
    work_directory = "/work_directory/"
    filtered_output_file = os.path.join(work_directory, "filtered_output_all.inf")
    filtered_coord_file = os.path.join(work_directory, "filtered_coord.out")
    path_to_output = os.path.join(work_directory, "filtered_ribo.fa")
    
    retrieve_ribo_coordinates(filtered_output_file, filtered_coord_file)

    parse_fasta(filtered_coord_file, path_to_fasta_dir, path_to_output)