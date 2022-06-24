
rule chromosome_notation:
    input:
        tumor = chr_notation_input,
        normal = normal_input,
    output:
        sample = os.path.join(config["output_dir"], "preprocessed_sample/{sample}.bam")
    params:
        normal_notation = normal_output
    log:
        os.path.join(config["log_dir"], "log/{sample}.log")
    shell:
        """
        notation=$(samtools view -H {input.tumor}|less|awk "NR == 2 {{print $2}}"|cut -f 2|cut -d ":" -f 2)
        if [ "$notation" = "Chr1" ] || [ "$notation" = "chr1" ]; then 
            samtools view -H {input.tumor}| sed  -e 's/SN:Chr\([0-9XY]*\)/SN:chr\\1/' -e 's/SN:ChrMT/SN:chrM/'|samtools reheader - {input.tumor} > {output}
        else
            samtools view -H {input.tumor} | sed  -e 's/SN:\([0-9XY]*\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' |samtools reheader - {input.tumor} > {output}
        fi


        notation_normal=$(samtools view -H {input.normal}|less|awk "NR == 2 {{print $2}}"|cut -f 2|cut -d ":" -f 2)
        if [ "$notation_normal" = "Chr1" ] || [ "$notation_normal" = "chr1" ]; then 
            samtools view -H {input.normal}| sed  -e 's/SN:Chr\([0-9XY]*\)/SN:chr\\1/' -e 's/SN:ChrMT/SN:chrM/'|samtools reheader - {input.normal} > {params.normal_notation}
        else
            samtools view -H {input.normal} | sed  -e 's/SN:\([0-9XY]*\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' |samtools reheader - {input.normal} > {params.normal_notation}
        fi

        2 > {log}

        """







