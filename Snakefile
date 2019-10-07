# Run RFMix v1.5.4 and optionally ADMIXTURE

# Variables and filepaths stored here
include: "scripts/snakemake_variables.py"

shell.executable("/usr/bin/bash")
shell.prefix("source ~/.bashrc; ")

rule all:
    input:
        expand(DATA_DIR + "anc{anc}.combined_global_anc_frac.txt", anc=ADMIX_POP_VAL)

rule create_admix_filter:
    input:
        DATA_DIR + ADMIX_METADATA
    output:
        temp(DATA_DIR + "anc{anc}_filter.txt")
    shell:
        """
        conda activate py36
        python {FILTER_SCRIPT} \
            --sample_data {input} \
            --pop_col {ADMIX_POP_COL} \
            --pop_val {wildcards.anc} \
            --out {output}
        conda deactivate
        """

rule filter_admix:
    input:
        dataset=DATA_DIR + ADMIX_DATA,
        filter=rules.create_admix_filter.output
    output:
        bcf=temp(DATA_DIR + "anc{anc}.filtered.bcf.gz"),
        idx=temp(DATA_DIR + "anc{anc}.filtered.bcf.gz.csi")
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
        map=DATA_DIR + "chr_map.txt"
    output:
        bcf=temp(DATA_DIR + "anc{anc}.filtered.anno.bcf.gz"),
        idx=temp(DATA_DIR + "anc{anc}.filtered.anno.bcf.gz.csi")
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
        bcf=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.filtered.anno.bcf.gz",
        idx=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.filtered.anno.bcf.gz.csi"
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
        bcf=temp(expand(DATA_DIR + "chr{{chr}}/anc{{anc}}/000{idx}.bcf", idx=[0,1])),
        idx=temp(expand(DATA_DIR + "chr{{chr}}/anc{{anc}}/000{idx}.bcf.csi", idx=[0,1])),
        readme=temp(DATA_DIR + "chr{chr}/anc{anc}/README.txt")
    params:
        out_dir=DATA_DIR + "chr{chr}/anc{anc}"
    shell:
        """
        mkdir -p {params.out_dir}
        conda activate bcftools-env
        bcftools isec -p {params.out_dir} -n=2 --collapse none -Ob {input.admix} {input.ref}
        conda deactivate
        """

rule merge:
    input:
        bcf=rules.intersect.output.bcf,
        idx=rules.intersect.output.idx
    output:
        temp(DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.merged.vcf.gz")
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
        keep=temp(DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.r2{r2}.prune.in"),
        exclude=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.r2{{r2}}.{ext}", ext=["log", "nosex", "prune.out"]))
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.r2{r2}"
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
        DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.merged.pruned.{r2}.vcf.gz"
    shell:
        """
        conda activate bcftools-env
        bcftools view -i 'ID=@{input.snps}' -Oz -o {output} {input.vcf}
        conda deactivate
        """

rule make_rfmix_input_acro:
    input:
        vcf=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.merged.pruned.{r2}.vcf.gz", r2=RFMIX_R2),
        admix_info=DATA_DIR + ADMIX_METADATA,
        ref_info=REF_DIR + REF_METADATA,
        genetic_map=MAP_DIR + "plink.chr{chr}.GRCh38.map"
    output:
        snp_locations=temp(DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.snp_locations.q.txt"),
        classes=temp(DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.classes.txt"),
        alleles=temp(DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.alleles.q.txt"),
        pop_map=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.pop_map.txt",
        pos_map=temp(DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.pos_map.q.txt")
    wildcard_constraints:
        chr="|".join(ACROCENTRIC)
    params:
        out_file="data/chr{chr}/chr{chr}.anc{anc}"
    shell:
        """
        conda activate py36
        CENT=0
        zcat {input.vcf} | python {RFMIX_INPUT_SCRIPT} \
            --admix {input.admix_info} \
            --admix_pop_col {ADMIX_POP_COL} \
            --admix_pop_val {wildcards.anc} \
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
        vcf=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.merged.pruned.{r2}.vcf.gz", r2=RFMIX_R2),
        admix_info=DATA_DIR + ADMIX_METADATA,
        ref_info=REF_DIR + REF_METADATA,
        genetic_map=MAP_DIR + "plink.chr{chr}.GRCh38.map",
        centromeres=MAP_DIR + "chrom_centromere.map"
    output:
        snp_locations=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.snp_locations.{arm}.txt", arm=ARMS)),
        classes=temp(DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.classes.txt"),
        alleles=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.alleles.{arm}.txt", arm=ARMS)),
        pop_map=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.pop_map.txt",
        pos_map=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.pos_map.{arm}.txt", arm=ARMS))
    wildcard_constraints:
        chr="|".join(NON_ACRO)
    params:
        out_file="data/chr{chr}/chr{chr}.anc{anc}"
    shell:
        """
        conda activate py36
        CENT=$(grep '^{wildcards.chr}[^0-9]' {input.centromeres} | cut -f3 -d$'\t')
        zcat {input.vcf} | python {RFMIX_INPUT_SCRIPT} \
            --admix {input.admix_info} \
            --admix_pop_col {ADMIX_POP_COL} \
            --admix_pop_val {wildcards.anc} \
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
        snp_locations=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.snp_locations.{arm}.txt",
        classes=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.classes.txt",
        alleles=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.alleles.{arm}.txt"
    output:
        temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.{{arm}}.{em}.Viterbi.txt", em=EM_ITER))
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.{arm}"
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
            "chr{{chr}}/chr{{chr}}.anc{{anc}}.{arm}.{{em}}.Viterbi.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.q.{em}.Viterbi.txt",
        snp_locations=lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.anc{{anc}}.snp_locations.{arm}.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.snp_locations.q.txt",
        pos_map=lambda wildcards: expand(DATA_DIR + \
            "chr{{chr}}/chr{{chr}}.anc{{anc}}.pos_map.{arm}.txt", arm=ARMS) \
            if wildcards.chr not in ACROCENTRIC else \
            DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.pos_map.q.txt",
    output:
        viterbi=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.em{em}.Viterbi.merged.txt",
        snp_locations=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.em{em}.snp_locations.merged.txt",
        pos_map=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.em{em}.pos_map.merged.txt"
    shell:
        """
        cat {input.viterbi} > {output.viterbi}
        cat {input.snp_locations} > {output.snp_locations}
        cat {input.pos_map} > {output.pos_map}
        """

rule collapse_ancestry:
    input:
        viterbi=rules.merge_rfmix_output.output.viterbi,
        pos_map=rules.merge_rfmix_output.output.pos_map,
        pop_map=rules.make_rfmix_input.output.pop_map
    output:
        tracts=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.em{em}.tracts.txt"
    shell:
        """
        conda activate py27
        python {RFMIX_OUTPUT_SCRIPT} \
            --rfmix {input.viterbi} \
            --snp_map {input.pos_map} \
            --pop_map {input.pop_map} \
            --pop_labels {REF_POP_VAL} \
            --out {output} \
            --chr {wildcards.chr}
        conda deactivate
        """

rule global_inference:
    input:
        tracts=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.em{em}.tracts.txt", em=EM_ITER)
    output:
        DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.lai_global.txt"
    shell:
        """
        conda activate py36
        python {GLOBAL_INF_SCRIPT} \
            --tracts {input.tracts} \
            --pops {REF_POP_VAL} \
            --out {output}
        conda deactivate
        """

rule plot_karyogram:
    input:
        tracts=expand(DATA_DIR + "chr{chr}/chr{chr}.anc{{anc}}.em{em}.tracts.txt", chr=CHROMS, em=EM_ITER),
        cent=MAP_DIR + "centromere.tsv"
    output:
        "plots/anc{anc}.{ind}.png"
    shell:
        """
        conda activate py36
        python {KARYOGRAM_SCRIPT} \
            --tracts {input.tracts} \
            --ind {wildcards.ind} \
            --centromeres {input.cent} \
            --pop_order {REF_POP_VAL} \
            --colors blue,red,black \
            --out {output}
        conda deactivate
        """

rule admixture_input:
    input:
        vcf=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.merged.pruned.{r2}.vcf.gz",
        map=rules.make_rfmix_input.output.pop_map
    output:
        bed=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.merged.pruned.{r2}.bed",
        misc=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.merged.pruned.{{r2}}.{ext}", ext=["bim", "fam"]),
        discard=temp(expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.merged.pruned.{{r2}}.{ext}", ext=["log", "nosex"])),
        pop=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.merged.pruned.{r2}.pop"
    params:
        out_file=DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.merged.pruned.{r2}"
    shell:
        """
        conda activate plink-env
        plink --vcf {input.vcf} --make-bed --keep-allele-order --out {params.out_file}
        conda deactivate
        tail -n 1 {input.map} | sed 's/\t/\n/g' | sed 's/ADMIX/-/g' > {output.pop}
        """

rule run_admixture:
    input:
        bed=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.merged.pruned.{r2}.bed", r2=ADMIX_R2),
        misc=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.merged.pruned.{r2}.{ext}", ext=["bim", "fam"], r2=ADMIX_R2),
        pop=expand(DATA_DIR + "chr{{chr}}/chr{{chr}}.anc{{anc}}.merged.pruned.{r2}.pop", r2=ADMIX_R2)
    output:
        expand("chr{{chr}}.anc{{anc}}.merged.pruned.{r2}.{n}.{ext}", ext=['P', 'Q'], n=NPOP, r2=ADMIX_R2)
    shell:
        """
        {ADMIXTURE} {input.bed} 2 --supervised
        """

rule combine_rfmix_admixture:
    input:
        rfmix=rules.global_inference.output,
        admix=expand("chr{{chr}}.anc{{anc}}.merged.pruned.{r2}.{n}.Q", n=NPOP, r2=ADMIX_R2),
        map=rules.make_rfmix_input.output.pop_map
    output:
        DATA_DIR + "chr{chr}/chr{chr}.anc{anc}.combined_global_anc_frac.txt"
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
        anc=expand(DATA_DIR + "chr{chr}/chr{chr}.anc{{anc}}.combined_global_anc_frac.txt", chr=CHROMS),
        map=MAP_DIR + "chrom_centromere.map"
    output:
        DATA_DIR + "anc{anc}.combined_global_anc_frac.txt"
    shell:
        """
        conda activate py36
        python {COMB_CHR_SCRIPT} \
            --anc {input.anc} \
            --map {input.map} \
            --out {output}
        conda deactivate
        """

# ### Below are temporary rules for handling Europeans
# rule create_eur_filter:
#     input:
#         DATA_DIR + ADMIX_METADATA
#     output:
#         temp(DATA_DIR + "eur_filter.txt")
#     shell:
#         """
#         conda activate py36
#         python {FILTER_SCRIPT} \
#             --sample_data {input} \
#             --pop_col {ADMIX_POP_COL} \
#             --pop_val 1 \
#             --out {output}
#         conda deactivate
#         """
#
# rule filter_eur:
#     input:
#         dataset=DATA_DIR + ADMIX_DATA,
#         filter=rules.create_eur_filter.output
#     output:
#         bcf=temp(DATA_DIR + ".eur.filtered.bcf.gz"),
#         idx=temp(DATA_DIR + ".eur.filtered.bcf.gz.csi")
#     shell:
#         """
#         conda activate bcftools-env
#         bcftools view -S {input.filter} --force-samples -Ob {input.dataset} | \
#             bcftools view -i 'MAF[0] > 0.05' -Ob | \
#             bcftools view --genotype ^miss --phased -Ob -o {output.bcf}
#         bcftools index -c {output.bcf}
#         conda deactivate
#         """
#

#
# rule split_eur:
#     input:
#         bcf=rules.annotate_eur.output.bcf,
#         idx=rules.annotate_eur.output.idx
#     output:
#         bcf=DATA_DIR + "chr{chr}/chr{chr}." + ".eur.filtered.anno.bcf.gz",
#         idx=DATA_DIR + "chr{chr}/chr{chr}." + ".eur.filtered.anno.bcf.gz.csi"
#     params:
#         out_dir=DATA_DIR + "chr{chr}"
#     shell:
#         """
#         conda activate bcftools-env
#         mkdir -p {params.out_dir}
#         bcftools view -r {wildcards.chr} -Ob -o {output.bcf} {input.bcf}
#         bcftools index -c {output.bcf}
#         conda deactivate
#         """
#
# rule intersect_eur:
#     input:
#         bcf=rules.split_eur.output.bcf,
#         idx=rules.split_eur.output.idx,
#         other=DATA_DIR + "chr{chr}/chr{chr}.merged.pruned.0.5.vcf.gz"
#     output:
#         bcf=temp(expand(DATA_DIR + "chr{{chr}}/000{idx}.eur.bcf", idx=[0])),
#         idx=temp(expand(DATA_DIR + "chr{{chr}}/000{idx}.eur.bcf.csi", idx=[0])),
#         readme=temp(DATA_DIR + "chr{chr}/README.txt")
#     params:
#         out_dir=DATA_DIR + "chr{chr}"
#     shell:
#         """
#         conda activate bcftools-env
#         bcftools isec -p {params.out_dir} -n=2 --collapse none -Ob {input.bcf} {input.other}
#         conda deactivate
#         """
