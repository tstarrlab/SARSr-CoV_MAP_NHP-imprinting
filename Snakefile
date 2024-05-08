"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import os
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_dag,
            make_summary,
            get_WH1_mut_antibody_escape,
            get_SARSr_wts_mut_bind_expr,
            get_SARSr_DMS_mut_bind_expr,
            get_SARS2_DMS_mut_bind_expr,
            save_pinned_env
            
            

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# Rules -----------------------------------------------------------------------

# making this summary is the target rule (in place of `all`) since it
# is first rule listed.
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
        env='environment_pinned.yml',
        SARSr_lib47_get_mut_bind_expr=config['SARSr_lib47_mut_bind_expr'],
        SARSr_lib40_get_mut_bind_expr=config['SARSr_lib40_mut_bind_expr'],
        SARS2_WH1_BA2_get_mut_bind_expr=config['SARS2_WH1_BA2_mut_bind_expr'],
        WH1_get_mut_antibody_escape=config['WH1_mut_antibody_escape'],
        variant_counts_file=config['variant_counts_file'],
        count_variants=nb_markdown('count_variants.ipynb'),
        compute_AUC='results/summary/compute_AUC.md',
        barcode_AUC=config['sera_delta_AUC_file'],
        collapse_bc_dms='results/summary/collapse_barcodes_SARSr-DMS.md',
        collapse_bc_dms_file=config['final_variant_scores_dms_file'],
        collapse_bc_wts='results/summary/collapse_barcodes_SARSr-wts.md',
        collapse_bc_wts_file=config['final_variant_scores_wts_file'],
        analyze_epitopes='results/summary/analyze_epitopes.md',
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the DAG of the computational workflow:
            ![{path(input.dag)}]({path(input.dag)})

            Here is the Markdown output of each Jupyter notebook in the
            workflow:



            1. Get prior Wuhan-Hu-1 antibody escape data from the [escape map aggregator repository](https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps/blob/main/processed_data/escape_data.csv), sarbecovirus homolog wildtypes expression measurements from [this repository](https://github.com/jbloomlab/SARSr-CoV_homolog_survey), sarbecovirus DMS library mutant expression measurements from [this repository](https://github.com/jbloomlab/SARSr-CoV-RBD_DMS), and SARS-CoV-2 WH1 and BA.2 DMS library mutant expression measurements from [this repository](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron). 
            
            2. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts_file)})
               giving counts of each barcoded variant in each condition.

            3. [Fit AUC to sera binding curves]({path(input.compute_AUC)}).
               Creates a [table]({path(input.barcode_AUC)})
               giving the AUC phenotype of each barcoded variant in each condition.
            
            4. Collapse internal replicate barcodes of each variant to final variant phenotypes for the wildtype sarbecovirus homologs pool. Analysis [here]({path(input.collapse_bc_wts)}) and final output file [here]({path(input.collapse_bc_wts_file)}).
            
            5. Collapse internal replicate barcodes of each variant to final variant phenotypes for the sarbecovirus DMS pools. Analysis [here]({path(input.collapse_bc_dms)}) and final output file [here]({path(input.collapse_bc_dms_file)}).
            
            6. [Analyze site-wise average binding values]({path(input.analyze_epitopes)}) to visualize the epitopes targeted by each vaccine/SARSr background combination. Outputs pdb files that can be used to visualize structural epitopes in PyMol.
                        
            """
            ).strip())

rule make_dag:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'dag.svg')
    shell:
        "snakemake --forceall --dag | dot -Tsvg > {output}"

rule save_pinned_env:
    input:
        dag=os.path.join(config['summary_dir'], 'dag.svg'),
    log:
    	"environment_pinned.yml"
    shell:
        """
        conda env export > {log}
        """

rule analyze_epitopes:
    input:
        config['final_variant_scores_dms_file'],
        config['WH1_mut_antibody_escape'],
    output:
        md='results/summary/analyze_epitopes.md',
        md_files=directory('results/summary/analyze_epitopes_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='analyze_epitopes.Rmd',
        md='analyze_epitopes.md',
        md_files='analyze_epitopes_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule collapse_bcs_SARSr_DMS:
    input:
        config['sera_delta_AUC_file'],
        config['SARSr_lib40_mut_bind_expr'],
        config['SARS2_WH1_BA2_mut_bind_expr'],
        config['RBD_annotation_file']
    output:
        config['final_variant_scores_dms_file'],
        md='results/summary/collapse_barcodes_SARSr-DMS.md',
        md_files=directory('results/summary/collapse_barcodes_SARSr-DMS_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='collapse_barcodes_SARSr-DMS.Rmd',
        md='collapse_barcodes_SARSr-DMS.md',
        md_files='collapse_barcodes_SARSr-DMS_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule collapse_bcs_SARSr_wts:
    input:
        config['sera_delta_AUC_file'],
        config['SARSr_lib47_mut_bind_expr']
    output:
        config['final_variant_scores_wts_file'],
        md='results/summary/collapse_barcodes_SARSr-wts.md',
        md_files=directory('results/summary/collapse_barcodes_SARSr-wts_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='collapse_barcodes_SARSr-wts.Rmd',
        md='collapse_barcodes_SARSr-wts.md',
        md_files='collapse_barcodes_SARSr-wts_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule calculate_bc_sera_AUC:
    input:
        config['barcode_runs'],
        config['codon_variant_table'],
        config['variant_counts_file']
    output:
        config['sera_delta_AUC_file'],
        md='results/summary/compute_AUC.md',
        md_files=directory('results/summary/compute_AUC_files')
    envmodules:
        'R/4.1.3'
    params:
        nb='compute_AUC.Rmd',
        md='compute_AUC.md',
        md_files='compute_AUC_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """


rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['barcode_runs'],
        config['codon_variant_table'],
    output:
        config['variant_counts_file'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule get_SARS2_DMS_mut_bind_expr:
    """Download SARS2 WH1 and BA2 DMS library ACE2-binding and expression scores from URL."""
    output:
        file=config['SARS2_WH1_BA2_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['SARS2_WH1_BA2_mut_bind_expr_url'], output.file)

rule get_SARSr_DMS_mut_bind_expr:
    """Download SARSr DMS library ACE2-binding and expression scores from URL."""
    output:
        file=config['SARSr_lib40_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['SARSr_lib40_mut_bind_expr_url'], output.file)
        
rule get_SARSr_wts_mut_bind_expr:
    """Download SARSr wts library ACE2-binding and expression scores from URL."""
    output:
        file=config['SARSr_lib47_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['SARSr_lib47_mut_bind_expr_url'], output.file)


rule get_WH1_mut_antibody_escape:
    """Download SARS-CoV-2 mutation antibody-escape data from URL."""
    output:
        file=config['WH1_mut_antibody_escape']
    run:
        urllib.request.urlretrieve(config['WH1_mut_antibody_escape_url'], output.file)
        
