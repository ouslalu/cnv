
rule create_file:
    input:
        
    output:
        os.path.join(config["input_dir"], "{sample}.txt"),
    log:
        os.path.join(config["log_dir"], "{sample}.log")
    shell:
        """
        echo {wildcards.sample} > {output}
        """
