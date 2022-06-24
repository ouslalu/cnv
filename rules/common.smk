import pandas as pd
import os


metadata = pd.read_csv(config["patient_metadata"])
metadata_df = metadata.set_index(["egan"], drop= False) 
tumor_df = metadata_df[(metadata_df["type"]== "a")|(metadata_df["type"]== "c")] 
normal_df = metadata_df[metadata_df["type"]== "b"]



def chr_notation_input(wildcards):
    patient_egan = metadata_df.loc[wildcards.sample]["egan"]
    egan_location = os.path.join(config["input_dir"],patient_egan+".bam")
    return egan_location


def tumor_sample(wildcards):
    tumor = tumor_df.loc[wildcards.sample]["egan"]
    patient_id = tumor_df.loc[wildcards.sample]["patient"]
    tumor_location = os.path.join(config["output_dir"],"preprocessed_sample", tumor +".bam")
    return tumor_location


def normal_input(wildcards):
    p_id = metadata_df.loc[wildcards.sample]["patient_id"]
    normal_egan = normal_df[normal_df["patient_id"]==p_id]["egan"][0]
    patient_normal_location = os.path.join(config["input_dir"], normal_egan+".bam")
    return patient_normal_location

def normal_output(wildcards):
    p_id = metadata_df.loc[wildcards.sample]["patient_id"]
    normal_egan = normal_df[normal_df["patient_id"]==p_id]["egan"][0]
    return os.path.join(config["output_dir"],normal_egan+ ".bam")



def rename_to_include_id(wildcards):
    patient_id = metadata_df.loc[wildcards.sample]["patient"]
    return patient_id


#os.path.join(config["output_dir"], "file_rename/{sample}.bam"