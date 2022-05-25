#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 06.06.2021
#@author: Nadezhda Sorokina
#@contact: nadya7n@bk.ru

import os
from tqdm import tqdm
from typing import Tuple


def parse_infernal_file(infernal_file: str) -> Tuple[str]:
    chunk_with_ribo = list()
    with open(infernal_file) as fh:
        infernal_file = fh.readlines()
        
        fact_absence_of_ribo = "No hits detected that satisfy reporting thresholds"
        
        if fact_absence_of_ribo not in infernal_file:
            for ind in range(len(infernal_file)):
                if infernal_file[ind].startswith(">>"):
                    header = infernal_file[ind].strip().split()
                    id_org = header[1]
                    desc_org = " ".join(header[2:])
                    character = infernal_file[ind+3].strip()
                    chunk_with_ribo.append((id_org, desc_org, character))
        return chunk_with_ribo
    
    
def parse_accession_and_query(infernal_file: str) -> Tuple[str]:
    with open(infernal_file) as fh:
        for line in fh:
            line = line.strip()
            if "Accession:" in line:
                accession = line.split()[1]
            if "Query:" in line:
                query = line.split()[1]
    return query, accession


def write_output_of_infernal(id_org: str, desc_org: str, query: str, accession: str, character: str, output_file: str) -> None:
    with open(output_file, "a+") as fw:
        fw.write(f">{id_org}:{query}:{accession}\t{desc_org}\n")
        fw.write(f"{character}\n")
        
        
def main(infernal_directory: str, output_directory: str) -> None:
    output_file = os.path.join(output_directory, "output.inf")
    for id_dir, obj in enumerate(tqdm(os.listdir(infernal_directory))):
        if obj.startswith("RF"):
            family_directory = os.path.join(infernal_directory, obj)
            for file in os.listdir(family_directory):
                file_name = os.path.join(infernal_directory, family_directory, file)
                query, accession = parse_accession_and_query(file_name)
                chunk_with_ribo = parse_infernal_file(file_name)
                if chunk_with_ribo:
                    for ribo in chunk_with_ribo:
                        id_org, desc_org, character = ribo
                        write_output_of_infernal(id_org, desc_org, query, accession, character, output_file)
        else:
            continue
    
    
if __name__ == "__main__":
    base_dir = "/base_dir/"
    infernal_directory = os.path.join(base_dir, "infernal_all_output")
    output_directory = os.path.join(base_dir, "parsed_infernal")
    
    main(infernal_directory, output_directory)
