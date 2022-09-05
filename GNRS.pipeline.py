import os


projectpath = config["projectpath"]
scriptpath = config["scriptpath"]
SAMPLES = config["SAMPLES"]




#usage: ~/anaconda3/envs/stainedglass/bin/snakemake -p -s GNRS.pipeline.py --configfile GNRS.pipeline.yaml --cores
rule all:
    input:
        projectpath + "/ONT/ONT.log",
        projectpath + "/ONT/ONT-cns.log",
        projectpath + "/ONT/ONT-map.log",
        projectpath + "/ONT/ONT-cns2.log",
        projectpath + "/ONT/index.log",
        projectpath + "/ONT/docker.log",
        projectpath + "/CLR/CLR-cns2.log",
        projectpath + "/CLR/lgs.log",
        projectpath + "/HiFi/HiFi.p_ctg.fa",
        projectpath + "/QUAST/quast.log",
        projectpath + "/extract_unmap/extract.log",
        projectpath + "/HiFi/per_fasta_for_AGE/seqkit.log",
        projectpath + "/age_realign/age.log",
        projectpath + "/age_realign/HiFi_two_end_anchor.txt",
        projectpath + "/age_realign/ONT_two_end_anchor.txt",
        projectpath + "/vcf_and_merge/vcf.log",
        projectpath + "/vcf_and_merge/representive/kalign.log",
        projectpath + "/vcf_and_merge/vcf_gz.log",
        projectpath + "/vg_graphaligner/vg_graph.log",
        projectpath + "/vg_giraffe/pack.log",
        projectpath + "/vg_graphaligner/741.ont.calls.Q0.vcf",
        projectpath + "/vg_graphaligner/JSC20.clr.calls.Q0.vcf",
        projectpath + "/vg_graphaligner/S288C.hifi.calls.Q0.vcf"




##########################################################
# ONT assembly
# mkdir output/ONT
rule wtdbg2_assembly:
    input:
        ont_long_read = config["ont_long_read"]
    output:
        assembly_log = projectpath + "/ONT/ONT.log",
        cns_log = projectpath + "/ONT/ONT-cns.log",
        map_log = projectpath + "/ONT/ONT-map.log",
        cns2_log = projectpath + "/ONT/ONT-cns2.log"
    params:
        wtdbg2 = config["wtdbg2"],
        wtpoa = config["wtpoa"],
        minimap2 = config["minimap2"],
        samtools = config["samtools"],
        assembly_path = projectpath + "/ONT/ONT",
        ctg_lay_gz = projectpath + "/ONT/ONT.ctg.lay.gz",
        ont_first = projectpath + "/ONT/ONT_first.fasta",
        first_bam = projectpath + "/ONT/ONT_first.sorted.bam",
        ont_two = projectpath + "/ONT/ONT_two.fasta"
    threads:
        1
    shell:
        "{params.wtdbg2} -i {input.ont_long_read}  -fo {params.assembly_path}  -t 64 -p 19 -AS 2 -s 0.05 -L 1000 -S 4 > {output.assembly_log} 2>&1"
        " && {params.wtpoa} -t 64  -i {params.ctg_lay_gz} -fo {params.ont_first} > {output.cns_log} 2>&1"
        " && {params.minimap2} -t 64 -ax map-ont -r2k {params.ont_first} {input.ont_long_read} | {params.samtools} sort -@ 10 > {params.first_bam} 2>{output.map_log}"
        " && {params.samtools} view --threads 10  -F0x900 {params.first_bam} | {params.wtpoa} -t 64 -d {params.ont_first} -i -  -fo {params.ont_two} > {output.cns2_log} 2>&1"
        

# ONT polishing
rule margin_polishing:
     input:
         ont_two = projectpath + "/ONT/ONT_two.fasta",
         ont_long_read = config["ont_long_read"]
     output:
         index_log = projectpath + "/ONT/index.log"
     params:
         minimap2 = config["minimap2"],
         samtools = config["samtools"],
         three_bam = projectpath + "/ONT/ONT_three.sorted.bam",
         three_filter_bam = projectpath + "/ONT/ONT_three.sorted.F104.bam"
     threads:
         1
     shell:
         "{params.minimap2} --MD -a -x map-ont -t 24 {input.ont_two} {input.ont_long_read} | {params.samtools} sort -@ 10 > {params.three_bam}"
         " && {params.samtools} view -@ 24 -hb -F 0x104 {params.three_bam} > {params.three_filter_bam}"
         " && {params.samtools} index {params.three_filter_bam} > {output.index_log}"


# ONT polishing
# cu02: need docker
rule dockering:
    input:
        ont_two = projectpath + "/ONT/ONT_two.fasta"
    output:
        docker_log = projectpath + "/ONT/docker.log"
    params:
        ont_path = projectpath + "/ONT",
        three_filter_bam = "ONT_three.sorted.F104.bam",
        ont_margin = "ONT_margin_polish",
        ont_two = "ONT_two.fasta"
    threads:
        1
    shell:
        "docker run -i --user=`id -u`:`id -g` --cpus='24' -v {params.ont_path}:/data kishwars/margin_polish:latest {params.three_filter_bam} {params.ont_two} "
        "/opt/MarginPolish/params/allParams.np.human.r94-g235.json -t 24 -o /data/{params.ont_margin} > {output.docker_log}"








##########################################################
# CLR assembly
# mkdir output/CLR
rule clr_assembly:
    input:
        clr_long_read = config["clr_long_read"]
    output:
        assembly_log = projectpath + "/CLR/CLR.log",
        cns_log = projectpath + "/CLR/CLR-cns.log",
        map_log = projectpath + "/CLR/CLR-map.log",
        cns2_log = projectpath + "/CLR/CLR-cns2.log"
    params:
        wtdbg2 = config["wtdbg2"],
        wtpoa = config["wtpoa"],
        minimap2 = config["minimap2"],
        samtools = config["samtools"],
        assembly_path = projectpath + "/CLR/CLR",
        ctg_lay_gz = projectpath + "/CLR/CLR.ctg.lay.gz",
        clr_first = projectpath + "/CLR/CLR_first.fasta",
        first_bam = projectpath + "/CLR/CLR_first.sorted.bam",
        clr_two = projectpath + "/CLR/CLR_two.fasta"
    threads:
        1
    shell:
        "{params.wtdbg2} -i {input.clr_long_read}  -fo {params.assembly_path}  -t 64 -p 19 -AS 2 -s 0.05 -L 1000 -S 4 > {output.assembly_log} 2>&1"
        " && {params.wtpoa} -t 64  -i {params.ctg_lay_gz} -fo {params.clr_first} > {output.cns_log} 2>&1"
        " && {params.minimap2} -t 64 -ax map-pb -r2k {params.clr_first} {input.clr_long_read} | {params.samtools} sort -@ 10 > {params.first_bam} 2>{output.map_log}"
        " && {params.samtools} view --threads 10  -F0x900 {params.first_bam} | {params.wtpoa} -t 64 -d {params.clr_first} -i -  -fo {params.clr_two} > {output.cns2_log} 2>&1"


# CLR polishing
rule nextpolishing:
    input:
        clr_two = projectpath + "/CLR/CLR_two.fasta",
        clr_long_read = config["clr_long_read"]
    output:
    	lgs_log = projectpath + "/CLR/lgs.log"
    params:
        minimap2 = config["minimap2"],
        samtools = config["samtools"],
        lgs_bam = projectpath + "/CLR/genome.lgs.bam",
        fofn = projectpath + "/CLR/pb.map.bam.fofn",
        python = config["python"],
        nextpolish = config["nextpolish"],
        lgs_fa = projectpath + "/CLR/genome.lgspolish.fa"
    threads:
        1
    shell:
        "{params.minimap2} -ax map-pb -t 24 {input.clr_two} {input.clr_long_read} | {params.samtools} sort - -m 2g --threads 24 -o {params.lgs_bam}"
        " && {params.samtools} index {params.lgs_bam} && ls {params.lgs_bam} > {params.fofn} && {params.python} {params.nextpolish} -g {input.clr_two} -l {params.fofn} -r clr -p 24 -sp -o {params.lgs_fa} > {output.lgs_log}"







##########################################################
# HiFi assembly
# mkdir output/HiFi
rule hifiasm_assembly:
	input:
	    hifi_long_read = config["hifi_long_read"]
	output:
	    hifi_fa = projectpath + "/HiFi/HiFi.p_ctg.fa"
	params:
	    hifiasm = config["hifiasm"],
	    hifi_out = projectpath + "/HiFi/HiFi",
	    hifi_gfa = projectpath + "/HiFi/HiFi.bp.p_ctg.gfa",
	    gfatools = config["gfatools"]
	threads:
	    1
	shell:
	    "{params.hifiasm} -z20 -o {params.hifi_out} -t 64 {input.hifi_long_read} && {params.gfatools} gfa2fa {params.hifi_gfa} > {output.hifi_fa}"
        










##########################################################
# QUAST get unmap region
# mkdir output/QUAST
# mkdir output/QUAST/CLR
# mkdir output/QUAST/ONT
# mkdir output/QUAST/HiFi
rule quast_evaluate:
	input:
	    fa_input_clr = projectpath + "/CLR/genome.lgspolish.fa",
	    fa_input_hifi = projectpath + "/HiFi/HiFi.p_ctg.fa",
	    fa_input_ont = projectpath + "/ONT/ONT_margin_polish.fa"
	output:
	    quast_log = projectpath + "/QUAST/quast.log"
	params:
	    quast = config["quast"],
	    dir_out_clr = projectpath + "/QUAST/CLR",
	    dir_out_hifi = projectpath + "/QUAST/HiFi",
	    dir_out_ont = projectpath + "/QUAST/ONT",
	    ref_genome = config["ref_genome"]
	threads:
	    1
	shell:
	    "{params.quast} --no-gc --no-plots --no-html --no-snps --min-contig 1000 -o {params.dir_out_clr} -r {params.ref_genome} -t 24 {input.fa_input_clr} > {output.quast_log} "
	    "&& {params.quast} --no-gc --no-plots --no-html --no-snps --min-contig 1000 -o {params.dir_out_hifi} -r {params.ref_genome} -t 24 {input.fa_input_hifi} > {output.quast_log} "
	    "&& {params.quast} --no-gc --no-plots --no-html --no-snps --min-contig 1000 -o {params.dir_out_ont} -r {params.ref_genome} -t 24 {input.fa_input_ont} > {output.quast_log}"


# mkdir output/extract_unmap
rule extract_unmap:
    input:
        unalign_info_clr = projectpath + "/QUAST/CLR/contigs_reports/contigs_report_genome-lgspolish.unaligned.info",
        unalign_info_hifi = projectpath + "/QUAST/HiFi/contigs_reports/contigs_report_HiFi-p_ctg.unaligned.info",
        unalign_info_ont = projectpath + "/QUAST/ONT/contigs_reports/contigs_report_ONT_margin_polish.unaligned.info",
        filtered_clr = projectpath + "/QUAST/CLR/contigs_reports/minimap_output/genome-lgspolish.coords.filtered",
        filtered_hifi = projectpath + "/QUAST/HiFi/contigs_reports/minimap_output/HiFi-p_ctg.coords.filtered",
        filtered_ont = projectpath + "/QUAST/ONT/contigs_reports/minimap_output/ONT_margin_polish.coords.filtered"
    output:
        extract_log = projectpath + "/extract_unmap/extract.log"
    params:
        python = config["python"],
        script1 = scriptpath + "/quast_unalign_coords_extract.py",
        nrs_clr = projectpath + "/extract_unmap/CLR_NRS_raw.txt",
        nrs_hifi = projectpath + "/extract_unmap/HiFi_NRS_raw.txt",
        nrs_ont = projectpath + "/extract_unmap/ONT_NRS_raw.txt"
    threads:
        1
    shell:
    	"{params.python} {params.script1} -quast_unalign_info {input.unalign_info_clr} -quast_coords {input.filtered_clr} -out_txt_NRS {params.nrs_clr} -sample_name CLR "
    	"&& {params.python} {params.script1} -quast_unalign_info {input.unalign_info_hifi} -quast_coords {input.filtered_hifi} -out_txt_NRS {params.nrs_hifi} -sample_name HiFi "
    	"&& {params.python} {params.script1} -quast_unalign_info {input.unalign_info_ont} -quast_coords {input.filtered_ont} -out_txt_NRS {params.nrs_ont} -sample_name ONT > {output.extract_log}"









##########################################################
# AGE Realign, get NRS positioning
# mkdir output/ONT/per_fasta_for_AGE
# mkdir output/CLR/per_fasta_for_AGE
# mkdir output/HiFi/per_fasta_for_AGE
rule get_per_fa:
    input:
        fa_input_hifi = projectpath + "/HiFi/HiFi.p_ctg.fa",
        fa_input_clr = projectpath + "/CLR/genome.lgspolish.fa",
	    fa_input_ont = projectpath + "/ONT/ONT_margin_polish.fa"
    output:
        seqkit_log = projectpath + "/HiFi/per_fasta_for_AGE/seqkit.log"
    params:
        seqkit = config["seqkit"],
        hifi_list = projectpath + "/HiFi/per_fasta_for_AGE/HiFi_name.list",
        hifi_list2 = projectpath + "/HiFi/per_fasta_for_AGE/HiFi_name.list2",
        clr_list = projectpath + "/CLR/per_fasta_for_AGE/CLR_name.list",
        clr_list2 = projectpath + "/CLR/per_fasta_for_AGE/CLR_name.list2",
        ont_list = projectpath + "/ONT/per_fasta_for_AGE/ONT_name.list",
        ont_list2 = projectpath + "/ONT/per_fasta_for_AGE/ONT_name.list2",
        out_fa_hifi = projectpath + "/HiFi/per_fasta_for_AGE/",
        out_fa_clr = projectpath + "/CLR/per_fasta_for_AGE/",
        out_fa_ont = projectpath + "/ONT/per_fasta_for_AGE/"
    threads:
        1
    shell:
        "{params.seqkit} seq --name {input.fa_input_hifi} -o {params.hifi_list} && awk '{{print $1}}' {params.hifi_list} > {params.hifi_list2} "
        "&& cat {params.hifi_list2} | xargs -t -P 20 -n1 -I {{}} {params.seqkit} grep -p {{}} {input.fa_input_hifi} -o {params.out_fa_hifi}HiFi_{{}}.fa > {output.seqkit_log} "
        "&& {params.seqkit} seq --name {input.fa_input_clr} -o {params.clr_list} && awk '{{print $1}}' {params.clr_list} > {params.clr_list2} "
        "&& cat {params.clr_list2} | xargs -t -P 20 -n1 -I {{}} {params.seqkit} grep -p {{}} {input.fa_input_clr} -o {params.out_fa_clr}CLR_{{}}.fa > {output.seqkit_log}  "
        "&& {params.seqkit} seq --name {input.fa_input_ont} -o {params.ont_list} && awk '{{print $1}}' {params.ont_list} > {params.ont_list2} "
        "&& cat {params.ont_list2} | xargs -t -P 20 -n1 -I {{}} {params.seqkit} grep -p {{}} {input.fa_input_ont} -o {params.out_fa_ont}ONT_{{}}.fa > {output.seqkit_log}  "


# mkdir output/age_realign
# mkdir output/age_realign/HiFi
# mkdir output/age_realign/ONT
# mkdir output/age_realign/CLR
# AGE: one of the commands exited with non-zero exit code
rule age_align:
    input:
        nrs_hifi = projectpath + "/extract_unmap/HiFi_NRS_raw.txt",
        nrs_clr = projectpath + "/extract_unmap/CLR_NRS_raw.txt",
        nrs_ont = projectpath + "/extract_unmap/ONT_NRS_raw.txt"
    output:
        age_log = projectpath + "/age_realign/age.log",
        hifi_txt = projectpath + "/age_realign/HiFi_two_end_anchor.txt",
        clr_txt = projectpath + "/age_realign/CLR_two_end_anchor.txt",
        ont_txt = projectpath + "/age_realign/ONT_two_end_anchor.txt"
    params:
        python = config["python"],
        script1 = scriptpath + "/batch_command_AGE.py",
        script2 = scriptpath + "/batch_parser_AGE_out.py",
        out_fa_hifi = projectpath + "/HiFi/per_fasta_for_AGE/",
        out_fa_clr = projectpath + "/CLR/per_fasta_for_AGE/",
        out_fa_ont = projectpath + "/ONT/per_fasta_for_AGE/",
        hifi_age = projectpath + "/age_realign/HiFi/HiFi_AGE_batch.sh",
        clr_age = projectpath + "/age_realign/CLR/CLR_AGE_batch.sh",
        ont_age = projectpath + "/age_realign/ONT/ONT_AGE_batch.sh"
    threads:
        1
    run:
        cmd1 = "cd output/age_realign/HiFi && %s %s -raw_txt %s -per_fasta_dir %s  -sample_name HiFi > %s && cat %s | xargs -n1 -t -P 24 -I {{}} bash {{}} > %s" % (params.python, params.script1, input.nrs_hifi, params.out_fa_hifi, params.hifi_age, params.hifi_age,  output.age_log)
        print(cmd1)
        os.system(cmd1)
        cmd2 = "cd output/age_realign/ && %s %s -age_batch_sh %s -sample_name HiFi > %s" % (params.python, params.script2, params.hifi_age, output.hifi_txt)
        print(cmd2)
        os.system(cmd2)
        cmd3 = "cd output/age_realign/CLR && %s %s -raw_txt %s -per_fasta_dir %s  -sample_name CLR > %s && cat %s | xargs -n1 -t -P 24 -I {{}} bash {{}} > %s" % (params.python, params.script1, input.nrs_clr, params.out_fa_clr, params.clr_age, params.clr_age,  output.age_log)
        print(cmd3)
        os.system(cmd3)
        cmd4 = "cd output/age_realign/ && %s %s -age_batch_sh %s -sample_name CLR > %s" % (params.python, params.script2, params.clr_age, output.clr_txt)
        print(cmd4)
        os.system(cmd4)
        cmd5 = "cd output/age_realign/ONT && %s %s -raw_txt %s -per_fasta_dir %s  -sample_name ONT > %s && cat %s | xargs -n1 -t -P 24 -I {{}} bash {{}} > %s" % (params.python, params.script1, input.nrs_ont, params.out_fa_ont, params.ont_age, params.ont_age,  output.age_log)
        print(cmd5)
        os.system(cmd5)
        cmd6 = "cd output/age_realign/ && %s %s -age_batch_sh %s -sample_name ONT > %s" % (params.python, params.script2, params.ont_age, output.ont_txt)
        print(cmd6)
        os.system(cmd6)








##########################################################
# VCF file of NRS
# mkdir output/vcf_and_merge
rule get_vcf:
    input:
        hifi_txt = projectpath + "/age_realign/HiFi_two_end_anchor.txt",
        clr_txt = projectpath + "/age_realign/CLR_two_end_anchor.txt",
        ont_txt = projectpath + "/age_realign/ONT_two_end_anchor.txt",
        fa_input_hifi = projectpath + "/HiFi/HiFi.p_ctg.fa",
        fa_input_clr = projectpath + "/CLR/genome.lgspolish.fa",
	    fa_input_ont = projectpath + "/ONT/ONT_margin_polish.fa"
    output:
        vcf_log = projectpath + "/vcf_and_merge/vcf.log"
    params:
        python = config["python"],
        script1 = scriptpath + "/txt_to_vcf.py",
        ref_genome = config["ref_genome"]
    threads:
        1
    shell:
        "cd output/vcf_and_merge && {params.python} {params.script1} -txt {input.hifi_txt} -assembly {input.fa_input_hifi} -sample_name HiFi -ref {params.ref_genome} > {output.vcf_log} "
        "&& {params.python} {params.script1} -txt {input.clr_txt} -assembly {input.fa_input_clr} -sample_name CLR -ref {params.ref_genome} > {output.vcf_log} "
        "&& {params.python} {params.script1} -txt {input.ont_txt} -assembly {input.fa_input_ont} -sample_name ONT -ref {params.ref_genome} > {output.vcf_log} "






##########################################################
# merge VCF file
# conda activate jasminesv_new
# mkdir output/vcf_and_merge/representive
# 3 is number of clique (SUPP >= 2 from jasmine)
rule jasmine_nrs:
    input:
        clr_vcf = projectpath + "/vcf_and_merge/CLR.final.vcf",
        ont_vcf = projectpath + "/vcf_and_merge/ONT.final.vcf",
        hifi_vcf = projectpath + "/vcf_and_merge/HiFi.final.vcf"
    output:
        kalign_log = projectpath + "/vcf_and_merge/representive/kalign.log",
        vcf_gz_log = projectpath + "/vcf_and_merge/vcf_gz.log"
    params:
        merge_vcf = projectpath + "/vcf_and_merge/jasmine.merge.vcf",
        python = config["python"],
        script1 = scriptpath + "/select_representive_by_kalign.py",
        dir_vcf = projectpath + "/vcf_and_merge/",
        parallel = config["parallel"],
        kalign = config["kalign"],
        script2 = scriptpath + "/kalign_afa_representive.py",
        script3 = scriptpath + "/get_final_representive_vcf.py",
        bgzip = config["bgzip"],
        tabix = config["tabix"]
    threads:
        1
    shell:
        "ls {input.clr_vcf} {input.ont_vcf} {input.hifi_vcf} > {params.dir_vcf}file_list.txt && jasmine file_list={params.dir_vcf}file_list.txt  out_file={params.merge_vcf} threads=10 --output_genotypes --ignore_strand --keep_var_ids "
        "&& cd output/vcf_and_merge/representive/ && {params.python} {params.script1} -vcf {params.merge_vcf} -dir {params.dir_vcf} -temp test "
        "&& (seq 1 3) | {params.parallel} -j 24  '{params.kalign} -i test_clique_{{}}.fa -o test_clique_{{}}.afa -f fasta' > {output.kalign_log} "
        "&& (seq 1 3) | {params.parallel} -j 24  '{params.python} {params.script2} -afa test_clique_{{}}.afa > test_clique_{{}}.score' "
        "&& cd ../ && {params.python} {params.script3} -dir_score representive/ -vcf {params.merge_vcf} -dir {params.dir_vcf} -out_name jasmine.merge.represent "
        "&& {params.bgzip} jasmine.merge.represent.final.vcf "
        "&& {params.tabix} -p vcf jasmine.merge.represent.final.vcf.gz > {output.vcf_gz_log}"










##########################################################
# vg graph genome
# mkdir output/vg_graphaligner 
rule vg_genome_graph:
    input:
        gz_vcf = projectpath + "/vcf_and_merge/jasmine.merge.represent.final.vcf.gz"
    output:
        vg_graph_log = projectpath + "/vg_graphaligner/vg_graph.log"
    params:
        vg = config["vg"],
        ref_genome = config["ref_genome"]
    threads:
        1
    shell:
        "cd output/vg_graphaligner && {params.vg} construct -a -f -p -S -r {params.ref_genome} -v {input.gz_vcf} -m 32 > x.vg "
        "&& {params.vg} index -t 24 -p -L -x x.xg x.vg "
        "&& {params.vg} gbwt -n 32 -g x.32.gg -o x.32.gbwt -x x.xg -P "
        "&& {params.vg} snarls -t 24 --include-trivial x.xg > x.trivial.snarls "
        "&& {params.vg} index -t 24 -j x.dist -x x.xg -s x.trivial.snarls "
        "&& {params.vg} minimizer -k 29 -w 11 -t 24 -i x.29.11.32.min -g x.32.gbwt -d x.dist x.xg > {output.vg_graph_log}"













##########################################################
# vg giraffe alignment
# mkdir output/vg_giraffe
rule vg_giraffe:
    input:
        pair_end_f1 = config["pair_end_f1"],
        pair_end_f2 = config["pair_end_f2"],
        gz_vcf = projectpath + "/vcf_and_merge/jasmine.merge.represent.final.vcf.gz",
        single_end = config["single_end"]
    output:
        pack_log = projectpath + "/vg_giraffe/pack.log"
    params:
        vg = config["vg"],
        x_xg = projectpath + "/vg_graphaligner/x.xg",
        x_32_gg = projectpath + "/vg_graphaligner/x.32.gg",
        x_min = projectpath + "/vg_graphaligner/x.29.11.32.min",
        x_dist = projectpath + "/vg_graphaligner/x.dist",
        x_gbwt = projectpath + "/vg_graphaligner/x.32.gbwt",
        gam_741 = projectpath + "/vg_giraffe/741.gam",
        pack_741 = projectpath + "/vg_giraffe/741.Q0.pack",
        snarls = projectpath + "/vg_graphaligner/x.trivial.snarls",
        vcf_741 = projectpath + "/vg_giraffe/741.NGS.calls.Q0.vcf",
        s288c_gam = projectpath + "/vg_giraffe/S288C.gam",
        s288c_pack = projectpath + "/vg_giraffe/S288C.Q0.pack",
        s288c_vcf = projectpath + "/vg_giraffe/S288C.NGS.calls.Q0.vcf"
    threads:
        1
    shell:
        "{params.vg} giraffe -x {params.x_xg} -g {params.x_32_gg} -m {params.x_min} -d {params.x_dist} -p -b default --rescue-algorithm dozeu -N ONT --gbwt-name {params.x_gbwt} -C 500 -o gam -f {input.pair_end_f1} -f {input.pair_end_f2} -t 24 > {params.gam_741} "
        "&& {params.vg} pack -x {params.x_xg} -g {params.gam_741} -Q 0 -t 24 -o {params.pack_741} > {output.pack_log} "
        "&& {params.vg} call {params.x_xg} -k {params.pack_741} -r {params.snarls} -t 24 -v {input.gz_vcf} -s 741 > {params.vcf_741} "
        "&& {params.vg} giraffe -x {params.x_xg} -g {params.x_32_gg} -m {params.x_min} -d {params.x_dist} -p -b default --rescue-algorithm dozeu -N S288C --gbwt-name {params.x_gbwt} -C 500 -o gam -f {input.single_end} -t 24 > {params.s288c_gam} "
        "&& {params.vg} pack -x {params.x_xg} -g {params.s288c_gam} -Q 0 -t 24 -o {params.s288c_pack} "
        "&& {params.vg} call {params.x_xg} -k {params.s288c_pack} -r {params.snarls} -t 24 -v {input.gz_vcf} -s S288C > {params.s288c_vcf}"






##########################################################
# graphaligner alignment
# ONT
rule graphaligner_align:
    input:
        gz_vcf = projectpath + "/vcf_and_merge/jasmine.merge.represent.final.vcf.gz",
        ont_long_read = config["ont_long_read"]
    output:
        vcf_741 = projectpath + "/vg_graphaligner/741.ont.calls.Q0.vcf"
    params:
        graphaligner = config["graphaligner"],
        vg = config["vg"],
        x_vg = projectpath + "/vg_graphaligner/x.vg",
        x_xg = projectpath + "/vg_graphaligner/x.xg",
        snarls = projectpath + "/vg_graphaligner/x.trivial.snarls",
        gam_741 = projectpath + "/vg_graphaligner/741.ont.gam",
        pack_741 = projectpath + "/vg_graphaligner/741.ont.Q0.pack"
    threads:
        1
    shell:
        "{params.graphaligner} -g {params.x_vg} -f {input.ont_long_read} -a {params.gam_741} -x vg -t 24 "
        "&& {params.vg} pack -x {params.x_xg} -g {params.gam_741} -Q 0 -t 24 -o {params.pack_741} "
        "&& {params.vg} call {params.x_xg} -k {params.pack_741} -r {params.snarls} -t 24 -v {input.gz_vcf} -s 741 > {output.vcf_741}"



# CLR
rule graphaligner_align2:
    input:
        gz_vcf = projectpath + "/vcf_and_merge/jasmine.merge.represent.final.vcf.gz",
        clr_long_read = config["clr_long_read"]
    output:
        vcf_jsc20 = projectpath + "/vg_graphaligner/JSC20.clr.calls.Q0.vcf"
    params:
        graphaligner = config["graphaligner"],
        vg = config["vg"],
        x_vg = projectpath + "/vg_graphaligner/x.vg",
        x_xg = projectpath + "/vg_graphaligner/x.xg",
        snarls = projectpath + "/vg_graphaligner/x.trivial.snarls",
        gam_jsc20 = projectpath + "/vg_graphaligner/JSC20.clr.gam",
        pack_jsc20 = projectpath + "/vg_graphaligner/JSC20.clr.Q0.pack"
    threads:
        1
    shell:
        "{params.graphaligner} -g {params.x_vg} -f {input.clr_long_read} -a {params.gam_jsc20} -x vg -t 24 "
        "&& {params.vg} pack -x {params.x_xg} -g {params.gam_jsc20} -Q 0 -t 24 -o {params.pack_jsc20} "
        "&& {params.vg} call {params.x_xg} -k {params.pack_jsc20} -r {params.snarls} -t 24 -v {input.gz_vcf} -s 741 > {output.vcf_jsc20}"




# HiFi
rule graphaligner_align3:
    input:
        gz_vcf = projectpath + "/vcf_and_merge/jasmine.merge.represent.final.vcf.gz",
        hifi_long_read = config["hifi_long_read"]
    output:
        vcf_s288c = projectpath + "/vg_graphaligner/S288C.hifi.calls.Q0.vcf"
    params:
        graphaligner = config["graphaligner"],
        vg = config["vg"],
        x_vg = projectpath + "/vg_graphaligner/x.vg",
        x_xg = projectpath + "/vg_graphaligner/x.xg",
        snarls = projectpath + "/vg_graphaligner/x.trivial.snarls",
        gam_s288c = projectpath + "/vg_graphaligner/S288C.hifi.gam",
        pack_s288c = projectpath + "/vg_graphaligner/S288C.hifi.Q0.pack"
    threads:
        1
    shell:
        "{params.graphaligner} -g {params.x_vg} -f {input.hifi_long_read} -a {params.gam_s288c} -x vg -t 24 "
        "&& {params.vg} pack -x {params.x_xg} -g {params.gam_s288c} -Q 0 -t 24 -o {params.pack_s288c} "
        "&& {params.vg} call {params.x_xg} -k {params.pack_s288c} -r {params.snarls} -t 24 -v {input.gz_vcf} -s 741 > {output.vcf_s288c}"



