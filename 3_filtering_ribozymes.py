#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 06.06.2021
#@author: Nadezhda Sorokina
#@contact: nadya7n@bk.ru

import os
from Bio import SeqIO
from typing import List


THRESHHOLD_LENGTHS = {"RF00008": 29, "RF00009": 152, "RF00010": 184, "RF00011": 183, "RF00028": 126, "RF00029": 47,
                      "RF00030": 133, "RF00094": 46, "RF00163": 23, "RF00173": 26, "RF00234": 86, "RF00373": 152,
                      "RF00621": 96, "RF00622": 39, "RF01577": 313, "RF01787": 41, "RF01788": 96, "RF01807": 100,
                      "RF01865": 71, "RF01998": 42, "RF01999": 59, "RF02001": 86, "RF02003": 69, "RF02004": 100,
                      "RF02005": 114, "RF02012": 74, "RF02275": 39, "RF02276": 34, "RF02277": 63, "RF02357": 104,
                      "RF02472": 60, "RF02678": 76, "RF02679": 35, "RF02681": 42, "RF02682": 37, "RF02684": 30,
                      "RF02685": 27, "RF02686": 76, "RF02688": 48, "RF03152": 39, "RF03154": 37, "RF03160": 34}


def write_filter_ribo(header: List[str], info: List[str], output_file: str) -> None:
    with open(output_file, "a+") as fw:
        header = "\t".join(header)
        info = "\t".join(info)
        fw.write(f">{header}\n{info}\n")
        
        
def filter_ribozymes_length(input_file: str, output_file: str) -> None:
    with open(input_file) as fh:
        parsed_file = fh.readlines()
        for ind in range(len(parsed_file)):
            if parsed_file[ind].startswith(">"):
                ribo_info = parsed_file[ind+1].strip().split()
                ribo_header = parsed_file[ind].strip().split("\t")
                
                # start_ribo index 9, stop_ribo index 10
                length_ribo = abs(int(ribo_info[10]) - int(ribo_info[9]))
                
                #id_org:name_of_family:id_family_ribo
                ribo_family = ribo_header[0].split(":")[2]
                
                if length_ribo >= TRASHHOLD_LENGTHS[ribo_family]:
                    write_filter_ribo(ribo_header, ribo_info, output_file)

                    
if __name__ == "__main__":
    base_path = "/base_path/"
    input_file = os.path.join(base_path, "output_all.inf")
    output_file = os.path.join(base_path, "filtered_output_all.inf")
    filter_ribozymes_length(input_file, output_file)