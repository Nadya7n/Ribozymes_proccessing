#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 06.06.2021
#@author: Nadezhda Sorokina
#@contact: nadya7n@bk.ru

import os 
from typing import List, Dict, 

def overlap_cds(path_input: str, path_output: str) -> None:
    RF_list, GENE_dict = list(), dict()
    with open(path_input) as fh:
        for line in fh:
            line = line.strip()
            if "RF" in line:
                RF_list.append(line)
            else:
                if not GENE_dict:
                    GENE_dict[1] = line 
                    if GENE_dict and RF_list:
                        GENE_dict[2] = line
                    chunk_1, chunk_2 = find_overlap(RF_list, GENE_dict) 
                    write_overlap_ribo(chunk_1, path_output)
                    write_overlap_ribo(chunk_2, path_output)
                    RF_list, GENE_dict = list(), dict()
                    else:
                        GENE_dict[1] = line

def find_overlap(RF_list: List[str], GENE_dict: Dict[int, str]) -> List[str]:
    chunk_1, chunk_2 = [GENE_dict[1]], [GENE_dict[2]]
    st_gene_1, end_gene_1 = start_end_coordinates(GENE_dict[1])
    st_gene_2, end_gene_2 = start_end_coordinates(GENE_dict[2])
    for element in RF_list:
        st_rf, end_rf = start_end_coordinates(element)
        chunk_1 = join_multiple_ribo(st_gene_1, end_gene_1, element, st_rf, chunk_1)
        chunk_2 = join_multiple_ribo(st_gene_2, end_gene_2, element, end_rf, chunk_2)
    return join_final_chunk(chunk_1), join_final_chunk(chunk_2)

def join_final_chunk(chunk: List[str]) -> str or None:
    if len(chunk) > 1:
        return "\t".join(chunk)

def join_multiple_ribo(s_gene: int, e_gene: int, ribo: str, coord_ribo: int, chunk: List[str]) -> List[str]:
    if coord_ribo in range(s_gene, e_gene + 1):
        chunk.append(ribo)
    return chunk

def start_end_coordinates(item: str) -> list:
    return [int(item) for item in item.split(":")[1:3]]

def write_overlap_ribo(chunk: str, path_to_output: str) -> None:
    if chunk:
        with open(path_to_output, "a+") as fw:
            fw.write(f"{chunk}\n")
            
            
def merge_with_gene_info(input_1: str, input_2: str, output: str) -> None:
    with open(input_1) as fh:
        dct_with_ribo_cds_overlap = dict()
        for line in fh:
            line = line.strip().split("\t")
            cds = line.pop(0)
            dct_with_ribo_cds_overlap[cds] = "\t".join(line) 
    with open(input_2) as fh_2:
        with open(output, "a+") as fw:
            for line in fh_2:
                line = line.strip().split("\t")

                if line[0] in dct_with_ribo_cds_overlap:
                    info_header = ":".join(line)
                    fw.write(f">{info_header}\n{dct_with_ribo_cds_overlap[line[0]]}\n")
            
            
if __name__ == "__main__":
    main_path = "/main_path/"
    path_input = os.path.join(main_path, "merge_cds_ribo.out")
    overlap_cds_ribo = os.path.join(main_path, "overlap_cds_ribo.out")
    info_gene_path = os.path.join(main_path, "info_cds_without_repeats.out")
    overlap_with_gene_names = os.path.join(main_path, "cds_names_ribo.out")
    overlap_cds(path_input, path_output)
    
    merge_with_gene_info(overlap_cds_ribo, info_gene_path, overlap_with_gene_names)
