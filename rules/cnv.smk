rule cnv:
    input:
        tumor=tumor_sample
    output:
        reference = os.path.join(config["output_dir"], "cnv/reference/{sample}.cnn")
    params:
        cns_output = directory(os.path.join(config["output_dir"], "cnv/{sample}")),
        normal = normal_output,
        genome = config["genome"],
        annotation = config["annotation"],
        average_size = config["average_target_size"],
        patient_id = rename_to_include_id
    threads:config["threads"]
    log:
        os.path.join(config["log_dir"], "cnv/{sample}.log")
    shell:
        """ 
        
        cnvkit.py batch {input.tumor} \
        -n {params.normal} -m wgs \
        -f {params.genome} \
        --target-avg-size {params.average_size}\
        --annotate {params.annotation} \
        --output-reference {output.reference} \
        --output-dir {params.cns_output}_{params.patient_id} \
        --diagram --scatter
        2 > {log}

        """

    

