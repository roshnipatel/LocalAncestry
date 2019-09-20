# Run RFMix v1.5.4 on MESA dataset (build 38)

# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

# Decide whether to validate local ancestry calls by comparing to ADMIXTURE
validate_admixture = True

if validate_admixture:
    rule all:
        input:
            expand("plots/{ind}.em0.png", ind=INDIV[:10])

rule create_admix_filter:
    input:
        DATA_DIR + ADMIX_METADATA
    output:
        temp(DATA_DIR + "admix_filter.txt")
    shell:
        """
        conda activate py36
        python {FILTER_SCRIPT} \
            --sample_data {input} \
            --pop_col {ADMIX_POP_COL} \
            --pop_val {ADMIX_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule filter_admix:
    input:
        dataset=DATA_DIR + ADMIX_DATA + ".vcf.gz",
        filter=rules.create_admix_filter.output
    output:
        bcf=temp(DATA_DIR + ADMIX_DATA + ".filter_admix.bcf.gz"),
        idx=temp(DATA_DIR + ADMIX_DATA + ".filter_admix.bcf.gz.csi")
    shell:
        """
        conda activate bcftools-env
        bcftools view -S {input.filter} --force-samples -Ob {input.dataset} | \
            bcftools view -i 'MAF[0] > 0.05' -Ob | \
            bcftools view --genotype ^miss --phased -Ob -o {output.bcf}
        bcftools index -c {output.bcf}
        conda deactivate
        """

# Renames MESA chromosomes for consistency with 1KG dataset.
rule annotate_admix:
    input:
        bcf=rules.filter_admix.output.bcf,
        idx=rules.filter_admix.output.idx,
        map=DATA_DIR + "chr_map.txt"
    output:
        bcf=temp(DATA_DIR + ADMIX_DATA + ".filter_admix.anno.bcf.gz"),
        idx=temp(DATA_DIR + ADMIX_DATA + ".filter_admix.anno.bcf.gz.csi")
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
        bcf=DATA_DIR + "chr{chr}/chr{chr}." + ADMIX_DATA + ".filter_admix.anno.bcf.gz",
        idx=DATA_DIR + "chr{chr}/chr{chr}." + ADMIX_DATA + ".filter_admix.anno.bcf.gz.csi",
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
        temp(DATA_DIR + "ref_filter.txt")
    shell:
        """
        conda activate py36
        python {FILTER_SCRIPT} \
            --sample_data {input} \
            --pop_col {REF_POP_COL} \
            --pop_val {REF_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule filter_ref:
    input:
        ref=REF_DIR + "chr{chr}/chr{chr}." + REF_DATA + ".vcf.gz",
        filter=rules.create_ref_filter.output
    output:
        bcf=DATA_DIR + "chr{chr}/chr{chr}." + REF_DATA + ".bcf.gz",
        idx=DATA_DIR + "chr{chr}/chr{chr}." + REF_DATA + ".bcf.gz.csi"
    params:
        out_dir=DATA_DIR + "chr{chr}"
    shell:
        """
        conda activate bcftools-env
        mkdir -p {params.out_dir}
        bcftools view -S {input.filter} --force-samples -Ob {input.ref} | \
            bcftools view --genotype ^miss --phased -Ob | \
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

rule merge:
    input:
        bcf=rules.intersect.output.bcf,
        idx=rules.intersect.output.idx
    output:
        temp(DATA_DIR + "chr{chr}/chr{chr}.merged.vcf.gz")
    shell:
        """
        conda activate bcftools-env
        bcftools merge -m none -Ob {input.bcf} | \
            bcftools view -i 'MAF[0] > .1' -Oz -o {output}
        conda deactivate
        """

# Input for this rule needs to be in VCF format (even though Plink can accept
# BCF format) because otherwise Plink throws a weird error that is, for whatever
# reason, not replicated in the VCF.
rule find_indep_snps:
    input:
        rules.merge.output
    output:
        keep=DATA_DIR + "chr{chr}/chr{chr}.{r2}.prune.in",
        exclude=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.{{r2}}.{ext}", ext=["log", "nosex", "prune.out"]))
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.{r2}"
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
        snps=DATA_DIR + "chr{chr}/chr{chr}.{r2}.prune.in"
    output:
        DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.{r2}.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'ID=@{input.snps}' -Oz -o {output} {input.vcf}
        conda deactivate
        """

rule make_rfmix_input_acro:
    input:
        vcf=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.0.5.vcf.gz",
        admix_info=DATA_DIR + ADMIX_METADATA,
        ref_info=REF_DIR + REF_METADATA,
        genetic_map=MAP_DIR + "plink.chr{chr}.GRCh38.map",
        centromeres=MAP_DIR + "centromeres.tsv"
    output:
        snp_locations=temp(DATA_DIR + "chr{chr}/chr{chr}.snp_locations.q.txt"),
        classes=temp(DATA_DIR + "chr{chr}/chr{chr}.classes.txt"),
        alleles=temp(DATA_DIR + "chr{chr}/chr{chr}.alleles.q.txt"),
        pop_map=DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt",
        pos_map=temp(DATA_DIR + "chr{chr}/chr{chr}.pos_map.q.txt")
    wildcard_constraints:
        chr="|".join(ACROCENTRIC)
    params:
        out_file="data/chr{chr}/chr{chr}"
    shell:
        """
        conda activate py36
        CENT=0
        zcat {input.vcf} | python {RFMIX_INPUT_SCRIPT} \
            --admix {input.admix_info} \
            --admix_pop_col {ADMIX_POP_COL} \
            --admix_pop_val {ADMIX_POP_VAL} \
            --ref {input.ref_info} \
            --ref_pop_col {REF_POP_COL} \
            --ref_pop_val {REF_POP_VAL} \
            --genetic_map {input.genetic_map} \
            --out {params.out_file} \
            --cent $CENT
        conda deactivate
        """

rule make_rfmix_input:
    input:
        vcf=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.0.5.vcf.gz",
        admix_info=DATA_DIR + ADMIX_METADATA,
        ref_info=REF_DIR + REF_METADATA,
        genetic_map=MAP_DIR + "plink.chr{chr}.GRCh38.map",
        centromeres=MAP_DIR + "centromeres.tsv"
    output:
        snp_locations=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.snp_locations.{arm}.txt", arm=ARMS)),
        classes=temp(DATA_DIR + "chr{chr}/chr{chr}.classes.txt"),
        alleles=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.alleles.{arm}.txt", arm=ARMS)),
        pop_map=DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt",
        pos_map=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.pos_map.{arm}.txt", arm=ARMS))
    wildcard_constraints:
        chr="|".join(NON_ACRO)
    params:
        out_file="data/chr{chr}/chr{chr}"
    shell:
        """
        conda activate py36
        CENT=$(grep 'CEN{wildcards.chr}[^0-9]' {input.centromeres} | cut -f3)
        zcat {input.vcf} | python {RFMIX_INPUT_SCRIPT} \
            --admix {input.admix_info} \
            --admix_pop_col {ADMIX_POP_COL} \
            --admix_pop_val {ADMIX_POP_VAL} \
            --ref {input.ref_info} \
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
        classes=rules.make_rfmix_input.output.classes,
        alleles=DATA_DIR + "chr{chr}/chr{chr}.alleles.{arm}.txt"
    output:
        temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.{{arm}}.{em}.Viterbi.txt", em=EM_ITER))
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
            "chr{{chr}}/chr{{chr}}.{arm}.{{em}}.Viterbi.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.q.{em}.Viterbi.txt",
        snp_locations=lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.snp_locations.{arm}.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.snp_locations.q.txt",
        pos_map=lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.pos_map.{arm}.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.pos_map.q.txt",
    output:
        viterbi=DATA_DIR + "chr{chr}/chr{chr}.{em}.Viterbi.merged.txt",
        snp_locations=DATA_DIR + "chr{chr}/chr{chr}.{em}.snp_locations.merged.txt",
        pos_map=DATA_DIR + "chr{chr}/chr{chr}.{em}.pos_map.merged.txt"
    shell:
        """
        cat {input.viterbi} > {output.viterbi}
        cat {input.snp_locations} > {output.snp_locations}
        cat {input.pos_map} > {output.pos_map}
        """

# rule pos_map:
#     input:
#         vcf=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.0.5.vcf.gz",
#         genetic_map=MAP_DIR + "plink.chr{chr}.GRCh38.map"
#     output:
#         DATA_DIR + "chr{chr}/chr{chr}.{em}.pos_map.merged.txt"
#     params:
#         out_file="data/chr{chr}/chr{chr}"
#     shell:
#         """
#         conda activate py36
#         zcat {input.vcf} | python scripts/pos_map.py \
#             --genetic_map {input.genetic_map} \
#             --out {params.out_file}
#         conda deactivate
#         """

# rule pop_map:
#     input:
#         admix_info=DATA_DIR + ADMIX_METADATA,
#         ref_info=REF_DIR + REF_METADATA,
#         samp=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.0.5.vcf.gz"
#     output:
#         DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt"
#     shell:
#         """
#         conda activate py36
#         python {POP_SCRIPT} \
#             --samples $(zgrep -m 1 '#CHROM' {input.samp} | cut -d$'\t' --complement -f1-9 | sed 's/\t/ /g') \
#             --admix {input.admix_info} \
#             --admix_pop_col {ADMIX_POP_COL} \
#             --admix_pop_val {ADMIX_POP_VAL} \
#             --ref {input.ref_info} \
#             --ref_pop_col {REF_POP_COL} \
#             --ref_pop_val {REF_POP_VAL} \
#             --out {output}
#         conda deactivate
#         """

rule collapse_ancestry:
    input:
        viterbi=rules.merge_rfmix_output.output.viterbi,
        pos_map=DATA_DIR + "chr{chr}/chr{chr}.0.pos_map.merged.txt",
        pop_map=DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt"
    output:
        bed=expand(DATA_DIR + "bed/em{{em}}.chr{{chr}}.{ind}.{hapl}.bed", hapl=HAPL, ind=INDIV),
    params:
        out_file=DATA_DIR + "bed/em{em}.chr{chr}",
        ind=" ".join(INDIV)
    shell:
        """
        conda activate py27
        python {RFMIX_OUTPUT_SCRIPT} \
            --rfmix {input.viterbi} \
            --snp_map {input.pos_map} \
            --ind {params.ind} \
            --pop_map {input.pop_map} \
            --pop_labels {REF_POP_VAL} \
            --out {params.out_file} \
            --chr {wildcards.chr}
        conda deactivate
        """

rule global_inference:
    input:
        bed=expand(DATA_DIR + "bed/em{em}.chr{{chr}}.{ind}.{hapl}.bed", ind=INDIV, hapl=HAPL, em=EM_ITER)
    output:
        DATA_DIR + "chr{chr}/chr{chr}.lai_global.txt"
    params:
        ind=" ".join(INDIV)
    shell:
        """
        conda activate py36
        python {GLOBAL_INF_SCRIPT} \
            --bed {input.bed} \
            --ind {params.ind} \
            --pops {REF_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule plot_karyogram:
    input:
        bed=expand(DATA_DIR + "bed/em{{em}}.chr{chr}.{{ind}}.{hapl}.bed", hapl=HAPL, chr=CHROMS),
        cent=MAP_DIR + "centromere.map"
    output:
        "plots/{ind}.em{em}.png"
    shell:
        """
        conda activate py36
        python {KARYOGRAM_SCRIPT} \
            --bed_path {input.bed} \
            --ind {wildcards.ind} \
            --centromeres {input.cent} \
            --pop_order {REF_POP_VAL} \
            --colors blue,red,black \
            --out {output}
        conda deactivate
        """

rule admixture_input:
    input:
        vcf=rules.prune.output,
        map=ancient(DATA_DIR + "chr{chr}/chr{chr}.pop_map.txt")
    output:
        plink=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.merged.pruned.{{r2}}.{ext}", ext=["bed", "bim", "fam"]),
        discard=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.merged.pruned.{{r2}}.{ext}", ext=["log", "nosex"])),
        pop=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.{r2}.pop"
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.{r2}"
    shell:
        """
        conda activate plink-env
        plink --vcf {input.vcf} --make-bed --keep-allele-order --out {params.out_file}
        conda deactivate
        tail -n 1 {input.map} | sed 's/\t/\n/g' | sed 's/ADMIX/-/g' > {output.pop}
        """

rule run_admixture:
    input:
        bed=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.0.1.bed",
        misc=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.merged.pruned.0.1.{ext}", ext=["pop", "bim", "fam"])
    output:
        expand(DATA_DIR + "chr{{chr}}.merged.pruned.0.1.{n}.{ext}", ext=['P', 'Q'], n=NPOP)
    shell:
        """
        {ADMIXTURE} {input.bed} 2 --supervised
        mv chr* {DATA_DIR}
        """

rule combine_rfmix_admixture:
    input:
        rfmix=rules.global_inference.output,
        admix=expand(DATA_DIR + "chr{{chr}}.merged.pruned.0.1.{n}.Q", n=NPOP),
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
