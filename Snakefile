# Run RFMix v1.5.4 and optionally ADMIXTURE. Used to generate Fig. 1a and S1.

# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

# Required to use conda environments
shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

rule all:
    input:
        DATA_DIR + "combined_global_anc_frac.txt",
        DATA_DIR + "individuals_to_exclude.txt",
        DATA_DIR + "european_tract_lengths.txt"

rule create_admix_filter:
    input:
        DATA_DIR + ADMIX_METADATA
    output:
        DATA_DIR + "admix_filter.txt"
    shell:
        """
        conda activate py36
        python {FILTER_SCRIPT} \
            --sample_data {input} \
            --id_col {ADMIX_ID_COL} \
            --pop_col {ADMIX_POP_COL} \
            --pop_val {ADMIX_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule filter_admix:
    input:
        dataset=DATA_DIR + ADMIX_DATA,
        filter=rules.create_admix_filter.output
    output:
        bcf=temp(DATA_DIR + "mesa.admix.filtered.bcf.gz"),
        idx=temp(DATA_DIR + "mesa.admix.filtered.bcf.gz.csi")
    shell:
        """
        conda activate bcftools-env
        bcftools view -S {input.filter} --force-samples -Ob {input.dataset} | \
            bcftools view -i 'MAF[0] > 0.05' -Ob | \
            bcftools view --genotype ^miss --phased -Ob -o {output.bcf}
        bcftools index -c {output.bcf}
        conda deactivate
        """

# Renames MESA dataset for consistency with 1KG dataset.
rule annotate_admix:
    input:
        bcf=rules.filter_admix.output.bcf,
        idx=rules.filter_admix.output.idx,
        map=MAP_DIR + "chr_map.txt"
    output:
        bcf=temp(DATA_DIR + "mesa.admix.filtered.anno.bcf.gz"),
        idx=temp(DATA_DIR + "mesa.admix.filtered.anno.bcf.gz.csi")
    shell:
        """
        conda activate bcftools-env
        bcftools annotate --rename-chrs {input.map} -Ob -o {output.bcf} \
            {input.bcf}
        bcftools index -c {output.bcf}
        conda deactivate
        """

rule split_admix:
    input:
        bcf=rules.annotate_admix.output.bcf,
        idx=rules.annotate_admix.output.idx
    output:
        bcf=DATA_DIR + "chr{chr}/chr{chr}.mesa.admix.filtered.anno.bcf.gz",
        idx=DATA_DIR + "chr{chr}/chr{chr}.mesa.admix.filtered.anno.bcf.gz.csi"
    params:
        out_dir=DATA_DIR + "chr{chr}"
    shell:
        """
        conda activate bcftools-env
        mkdir -p {params.out_dir}
        bcftools view -r {wildcards.chr} -Ob -o {output.bcf} {input.bcf}
        bcftools index -c {output.bcf}
        conda deactivate
        """

rule create_ref_filter:
    input:
        REF_DIR + REF_METADATA
    output:
        DATA_DIR + "ref_filter.txt"
    shell:
        """
        conda activate py36
        python {FILTER_SCRIPT} \
            --sample_data {input} \
            --id_col {REF_ID_COL} \
            --pop_col {REF_POP_COL} \
            --pop_val {REF_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule filter_ref:
    input:
        ref=REF_DIR + "chr{chr}/chr{chr}." + REF_DATA,
        filter=rules.create_ref_filter.output
    output:
        bcf=DATA_DIR + "chr{chr}/chr{chr}.reference.bcf.gz",
        idx=DATA_DIR + "chr{chr}/chr{chr}.reference.bcf.gz.csi"
    params:
        out_dir=DATA_DIR + "chr{chr}"
    shell:
        """
        conda activate bcftools-env
        mkdir -p {params.out_dir}
        bcftools view -S {input.filter} --force-samples -Ou {input.ref} | \
            bcftools view --genotype ^miss --phased -Ou | \
            bcftools view -i 'MAF[0] > .05' -Ob -o {output.bcf}
        bcftools index -c {output.bcf}
        conda deactivate
        """

rule intersect:
    input:
        admix=rules.split_admix.output.bcf,
        admix_idx=rules.split_admix.output.idx,
        ref=rules.filter_ref.output.bcf,
        ref_idx=rules.filter_ref.output.idx
    output:
        bcf=temp(expand(DATA_DIR + "chr{{chr}}/000{idx}.bcf", idx=[0,1])),
        idx=temp(expand(DATA_DIR + "chr{{chr}}/000{idx}.bcf.csi", idx=[0,1])),
        readme=temp(DATA_DIR + "chr{chr}/README.txt")
    params:
        out_dir=DATA_DIR + "chr{chr}"
    shell:
        """
        conda activate bcftools-env
        bcftools isec -p {params.out_dir} -n=2 --collapse none -Ob {input.admix} {input.ref}
        conda deactivate
        """

# Explicitly filters for biallelic SNPs.
rule merge:
    input:
        bcf=rules.intersect.output.bcf,
        idx=rules.intersect.output.idx
    output:
        temp(DATA_DIR + "chr{chr}/chr{chr}.merged.maf{maf,0.[0-9]+}.vcf.gz")
    shell:
        """
        conda activate bcftools-env
        bcftools merge -m none -Ou {input.bcf} | \
            bcftools view -i 'MAF[0] > {wildcards.maf}' -m2 -M2 -v snps -Oz -o {output}
        conda deactivate
        """

# Input for this rule needs to be in VCF format (even though Plink can accept
# BCF format) because otherwise Plink throws a weird error that is, for whatever
# reason, not replicated in the VCF.
rule find_indep_snps:
    input:
        rules.merge.output
    output:
        keep=temp(DATA_DIR + "chr{chr}/chr{chr}.maf{maf}.r2{r2}.prune.in"),
        exclude=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.maf{{maf}}.r2{{r2}}.{ext}", 
                     ext=["log", "nosex", "prune.out"]))
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.maf{maf}.r2{r2}"
    shell:
        """
        conda activate plink-env
        plink --vcf {input} --vcf-half-call m --allow-no-sex \
            --keep-allele-order --indep-pairwise 50 5 {wildcards.r2} --out {params.out_file}
        conda deactivate
        """

rule prune:
    input:
        vcf=rules.merge.output,
        snps=rules.find_indep_snps.output.keep
    output:
        DATA_DIR + "chr{chr}/chr{chr}.merged.maf{maf,0.[0-9]+}.pruned.r2{r2, 0.[0-9]+}.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'ID=@{input.snps}' -Oz -o {output} {input.vcf}
        conda deactivate
        """

rule make_rfmix_input_acro:
    input:
        vcf=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.merged.maf{maf}.pruned.r2{r2}.vcf.gz", 
                   r2=RFMIX_R2, maf=RFMIX_MAF),
        admix_info=DATA_DIR + ADMIX_METADATA,
        ref_info=REF_DIR + REF_METADATA,
        genetic_map=MAP_DIR + "plink.chr{chr}.GRCh38.map"
    output:
        snp_locations=DATA_DIR + "chr{chr}/chr{chr}.snp_locations.q.txt",
        classes=DATA_DIR + "chr{chr}/chr{chr}.classes.txt",
        alleles=DATA_DIR + "chr{chr}/chr{chr}.alleles.q.txt",
        pop_map=DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt",
        pos_map=DATA_DIR + "chr{chr}/chr{chr}.pos_map.q.txt"
    wildcard_constraints:
        chr="|".join(ACROCENTRIC)
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}"
    shell:
        """
        conda activate py36
        CENT=0
        zcat {input.vcf} | python {RFMIX_INPUT_SCRIPT} \
            --admix {input.admix_info} \
            --admix_id_col {ADMIX_ID_COL} \
            --admix_pop_col {ADMIX_POP_COL} \
            --admix_pop_val {ADMIX_POP_VAL} \
            --ref {input.ref_info} \
            --ref_id_col {REF_ID_COL} \
            --ref_pop_col {REF_POP_COL} \
            --ref_pop_val {REF_POP_VAL} \
            --genetic_map {input.genetic_map} \
            --out {params.out_file} \
            --cent $CENT
        conda deactivate
        """

rule make_rfmix_input:
    input:
        vcf=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.merged.maf{maf}.pruned.r2{r2}.vcf.gz", 
                   r2=RFMIX_R2, maf=RFMIX_MAF),
        admix_info=DATA_DIR + ADMIX_METADATA,
        ref_info=REF_DIR + REF_METADATA,
        genetic_map=MAP_DIR + "plink.chr{chr}.GRCh38.map",
        centromeres=MAP_DIR + "chrom_centromere.map"
    output:
        snp_locations=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.snp_locations.{arm}.txt", arm=ARMS),
        classes=DATA_DIR + "chr{chr}/chr{chr}.classes.txt",
        alleles=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.alleles.{arm}.txt", arm=ARMS),
        pop_map=DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt",
        pos_map=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.pos_map.{arm}.txt", arm=ARMS)
    wildcard_constraints:
        chr="|".join(NON_ACRO)
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}"
    shell:
        """
        conda activate py36
        CENT=$(grep '^{wildcards.chr}[^0-9]' {input.centromeres} | cut -f3 -d$'\t')
        zcat {input.vcf} | python {RFMIX_INPUT_SCRIPT} \
            --admix {input.admix_info} \
            --admix_id_col {ADMIX_ID_COL} \
            --admix_pop_col {ADMIX_POP_COL} \
            --admix_pop_val {ADMIX_POP_VAL} \
            --ref {input.ref_info} \
            --ref_id_col {REF_ID_COL} \
            --ref_pop_col {REF_POP_COL} \
            --ref_pop_val {REF_POP_VAL} \
            --genetic_map {input.genetic_map} \
            --out {params.out_file} \
            --cent $CENT
        conda deactivate
        """

rule run_rfmix:
    input:
        snp_locations=DATA_DIR + "chr{chr}/chr{chr}.snp_locations.{arm}.txt",
        classes=DATA_DIR + "chr{chr}/chr{chr}.classes.txt",
        alleles=DATA_DIR + "chr{chr}/chr{chr}.alleles.{arm}.txt"
    output:
        DATA_DIR + "chr{chr}/chr{chr}.{arm}.0.Viterbi.txt",
        DATA_DIR + "chr{chr}/chr{chr}.{arm}.allelesRephased0.txt"
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.{arm}"
    shell:
        """
        conda activate py27
        python {RFMIX} PopPhased {input.alleles} \
            {input.classes} {input.snp_locations} \
            -e {EM_ITER} \
            -o {params.out_file} \
            --generations 8
        conda deactivate
        """

rule merge_rfmix_output:
    input:
        viterbi=lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.{arm}.0.Viterbi.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.q.0.Viterbi.txt",
        snp_locations=lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.snp_locations.{arm}.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.snp_locations.q.txt",
        pos_map=lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.pos_map.{arm}.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.pos_map.q.txt",
    output:
        viterbi=DATA_DIR + "chr{chr}/chr{chr}.Viterbi.merged.txt",
        snp_locations=DATA_DIR + "chr{chr}/chr{chr}.snp_locations.merged.txt",
        pos_map=DATA_DIR + "chr{chr}/chr{chr}.pos_map.merged.txt"
    shell:
        """
        cat {input.viterbi} > {output.viterbi}
        cat {input.snp_locations} > {output.snp_locations}
        paste {input.pos_map} > {output.pos_map}
        """

checkpoint collapse_ancestry:
    input:
        viterbi=rules.merge_rfmix_output.output.viterbi,
        pos_map=rules.merge_rfmix_output.output.pos_map,
        pop_map=rules.make_rfmix_input.output.pop_map
    output:
        completion_file=DATA_DIR + "bed/chr{chr}/job_complete.txt"
    params:
        out_dir=DATA_DIR + "bed/chr{chr}/",
    shell:
        """
        mkdir -p {params.out_dir}
        conda activate py27
        python {RFMIX_OUTPUT_SCRIPT} \
            --rfmix {input.viterbi} \
            --snp_map {input.pos_map} \
            --pop_map {input.pop_map} \
            --pop_labels {REF_POP_VAL} \
            --out {params.out_dir} \
            --chrom {wildcards.chr}
        conda deactivate
        """

def aggregate_bed_files(wildcards):
    checkpoint_output = checkpoints.collapse_ancestry.get(**wildcards).output[0]
    return(expand(DATA_DIR + "bed/chr{chr}/{ind}.{hapl}.bed", chr=wildcards.chr, \
                  ind=glob_wildcards(os.path.join(checkpoint_output, "{ind}.{hapl}.bed")).ind, \
                  hapl=glob_wildcards(os.path.join(checkpoint_output, "{ind}.{hapl}.bed")).hapl))
                  

rule global_inference:
    input:
        DATA_DIR + "bed/chr{chr}/job_complete.txt"
    output:
        DATA_DIR + "chr{chr}/chr{chr}.lai_global.txt"
    params:
        bed=DATA_DIR + "bed/chr{chr}/"
    shell:
        """
        conda activate py36
        python {GLOBAL_INF_SCRIPT} \
            --bed {params.bed} \
            --pops {REF_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule plot_karyogram:
    input:
        bed=expand(DATA_DIR + "bed/chr{chr}/{{ind}}.{hapl}.bed", hapl=HAPL, chr=CHROMS),
        cent=MAP_DIR + "centromere.tsv"
    output:
        "data/{ind}.png"
    shell:
        """
        conda activate py36
        python {KARYOGRAM_SCRIPT} \
            --bed_path {input.bed} \
            --ind {wildcards.ind} \
            --centromeres {input.cent} \
            --pop_order {REF_POP_VAL} \
            --colors {HEX_COLORS} \
            --out {output}
        conda deactivate
        """

rule admixture_pop_input:
    input:
        map=DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt"
    output:
        pop=DATA_DIR + "chr{chr}/chr{chr}.admixture.pop"
    shell:
        """
        tr '\t' '\n' < <(head -n 2 {input.map} | tail -n 1) | sed 's/ADMIX/-/' > {output.pop}
        """

rule admixture_snp_input:
    input:
        vcf=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.merged.maf{maf}.pruned.r2{r2}.vcf.gz", 
                   r2=ADMIX_R2, maf=ADMIX_MAF)
    output:
        bed=temp(DATA_DIR + "chr{chr}/chr{chr}.admixture.bed"),
        misc=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.admixture.{ext}", ext=["bim", "fam"])),
        discard=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.admixture.{ext}", ext=["log", "nosex"]))
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.admixture"
    shell:
        """
        conda activate plink-env
        plink --vcf {input.vcf} --make-bed --keep-allele-order --out {params.out_file}
        conda deactivate
        """

rule run_admixture:
    input:
        bed=DATA_DIR + "chr{chr}/chr{chr}.admixture.bed",
        misc=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.admixture.{ext}", ext=["bim", "fam"]),
        pop=DATA_DIR + "chr{chr}/chr{chr}.admixture.pop"
    output:
        expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.admixture.{n}.{ext}", ext=['P', 'Q'], n=NPOP)
    shell:
        """
        {ADMIXTURE} {input.bed} 2 --supervised
        mv chr{wildcards.chr}* {DATA_DIR}chr{wildcards.chr}/
        """

rule combine_rfmix_admixture:
    input:
        rfmix=rules.global_inference.output,
        admix=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.admixture.{n}.Q", 
                     n=NPOP, r2=ADMIX_R2),
        map=DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt"
    output:
        DATA_DIR + "chr{chr}/chr{chr}.combined_global_anc_frac.txt"
    shell:
        """
        conda activate py36
        python {COMB_SCRIPT} \
            --rfmix {input.rfmix} \
            --admix {input.admix} \
            --map {input.map} \
            --pops {REF_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule combine_chrs:
    input:
        anc=expand(rules.combine_rfmix_admixture.output, chr=CHROMS),
        map=MAP_DIR + "chrom_centromere.map"
    output:
        frac=DATA_DIR + "combined_global_anc_frac.txt",
        idv=DATA_DIR + "individuals_to_exclude.txt"
    shell:
        """
        conda activate py36
        python {COMB_CHR_SCRIPT} \
            --anc {input.anc} \
            --map {input.map} \
            --frac_out {output.frac} \
            --idv_out {output.idv} 
        conda deactivate
        """

rule merge_phasing:
    input:
        lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.{arm}.allelesRephased0.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.q.allelesRephased0.txt",
    output:
        DATA_DIR + "chr{chr}/chr{chr}.allelesRephased0.txt"
    shell:
        """
        cat {input} > {output}
        """

rule check_phasing:
    input:
        vcf=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.merged.pruned.{r2}.vcf.gz", 
                   r2=RFMIX_R2),
        phasing=DATA_DIR + "chr{chr}/chr{chr}.allelesRephased0.txt"
    output:
        DATA_DIR + "chr{chr}/chr{chr}.phasing_check.txt"
    shell:
        """
        conda activate py36
        python {PHASING_SCRIPT} \
            --vcf {input.vcf} \
            --rfmix_phased {input.phasing} \
            --out {output}
        conda deactivate
        """

rule compute_tract_lengths:
    input:
        expand(DATA_DIR + "bed/chr{chr}/job_complete.txt", chr=CHROMS)
    output:
        tract=DATA_DIR + "european_tract_lengths.txt"
    params:
        bed=DATA_DIR + "bed/"
    shell:
        """
        conda activate py36
        python {TRACT_SCRIPT} \
            --bed {params.bed} \
            --tract_out {output.tract} \
            --idv_out {output.idv}
        conda deactivate
        """
