package Genome::Model::ClinSeq;

use strict;
use warnings;
use Genome;
use List::MoreUtils;

# these are used below, and are also used in the documentation on commands in the tree
# to provide the most useful examples possible
our $DEFAULT_CANCER_ANNOTATION_DB_ID    = 'tgi/cancer-annotation/human/build37-20150205.1';
our $DEFAULT_MISC_ANNOTATION_DB_ID      = 'tgi/misc-annotation/human/build37-20130113.1';
our $DEFAULT_COSMIC_ANNOTATION_DB_ID    = 'cosmic/65.3';

class Genome::Model::ClinSeq {
    is => ['Genome::Model',
           'Genome::Model::ClinSeq::Util',],
    has_optional_input => [
        #Should there be a "somtic" base class for SomVal and SomVar?
        wgs_model               => { is => 'Genome::Model', doc => 'somatic model for wgs data' },
        exome_model             => { is => 'Genome::Model', doc => 'somatic model for exome data' },
        
        tumor_rnaseq_model      => { is => 'Genome::Model::RnaSeq', doc => 'rnaseq model for tumor rna-seq data' },
        normal_rnaseq_model     => { is => 'Genome::Model::RnaSeq', doc => 'rnaseq model for normal rna-seq data' },
        de_model                => { is => 'Genome::Model::DifferentialExpression', doc => 'differential-expression for tumor vs normal rna-seq data' },

        cancer_annotation_db    => { is => 'Genome::Db::Tgi::CancerAnnotation', default_value => $DEFAULT_CANCER_ANNOTATION_DB_ID },
        misc_annotation_db      => { is => 'Genome::Db::Tgi::MiscAnnotation', default_value => $DEFAULT_MISC_ANNOTATION_DB_ID },
        cosmic_annotation_db    => { is => 'Genome::Db::Cosmic', default_value => $DEFAULT_COSMIC_ANNOTATION_DB_ID },

        force                   => { is => 'Boolean', doc => 'skip sanity checks on input models' },

        #processing_profile      => { is => 'Genome::ProcessingProfile::ClinSeq', id_by => 'processing_profile_id', default_value => { } },
    ],
    #Processing profile parameters
    #If you add/modify params here make sure to add/modify params in
    # Genome/Test/Factory/ProcessingProfile/ClinSeq.pm
    has_param => [
        bam_readcount_version => { is => 'Text', doc => 'The bam readcount version to use during clonality analysis' },
        sireport_min_tumor_vaf => { is => 'Number', doc => 'Variants with a tumor VAF less than this (in any tumor sample) will be filtered out.'},
        sireport_max_normal_vaf => { is => 'Number', doc => 'Variants with a normal VAF greater than this (in any normal sample) will be filtered out.'},
        sireport_min_coverage => { is => 'Number', doc => 'Variants with coverage less than this (in any sample) will be filtered out.'},
        sireport_min_mq_bq => { is => 'Text', doc => 'Comma separated increasing-list of minimum mapping qualities for bam-readcounting.'},
        exome_cnv => { is => 'Boolean', doc => 'Should we run the Exome-cnv step or not.'},
    ],
    has_optional_metric => [
        common_name         => { is => 'Text', doc => 'the name chosen for the root directory in the build' },
    ],
    has_calculated => [
        expected_common_name => {
            is => 'Text',
            calculate_from => [qw /wgs_model exome_model tumor_rnaseq_model normal_rnaseq_model/],
            calculate => q|
              my ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
              if ($wgs_model) {
                  my $wgs_model_subject = $wgs_model->subject;
                  my $wgs_model_individual;
                  if ($wgs_model_subject->class eq 'Genome::Individual') {
                      $wgs_model_individual = $wgs_model_subject;
                  } else {
                      $wgs_model_individual = $wgs_model->subject->individual;
                  }
                  $wgs_common_name = $wgs_model_individual->common_name;
                  $wgs_name = $wgs_model_individual->name;
              }
              if ($exome_model) {
                  my $exome_model_subject = $exome_model->subject;
                  my $exome_model_individual;
                  if ($exome_model_subject->class eq 'Genome::Individual') {
                      $exome_model_individual = $exome_model_subject;
                  } else {
                      $exome_model_individual = $exome_model->subject->individual;
                  }
                  $exome_common_name = $exome_model_individual->common_name;
                  $exome_name = $exome_model_individual->name;
              }
              if ($tumor_rnaseq_model) {
                  $tumor_rnaseq_common_name = $tumor_rnaseq_model->subject->individual->common_name;
                  $tumor_rnaseq_name = $tumor_rnaseq_model->subject->individual->name;
              }
              if ($normal_rnaseq_model) {
                  $normal_rnaseq_common_name = $normal_rnaseq_model->subject->individual->common_name;
                  $normal_rnaseq_name = $normal_rnaseq_model->subject->individual->name;
              }
              my @names = ($wgs_common_name, $exome_common_name, $tumor_rnaseq_common_name, $normal_rnaseq_common_name, $wgs_name, $exome_name, $tumor_rnaseq_name, $normal_rnaseq_name);
              my $final_name = "UnknownName";
              foreach my $name (@names){
                if ($name){
                  $final_name = $name;
                  last();
                }
              }
              return $final_name;
            |
        },
    ],
    has => [
        individual_common_name => {
            via => '__self__',
            to => 'expected_common_name',
            doc => 'alias for benefit of Solr indexing',
        },
    ],
    doc => 'clinical and discovery sequencing data analysis and convergence of RNASeq, WGS and exome capture data',
};

sub define_by { 'Genome::Model::Command::Define::BaseMinimal' }

sub _help_synopsis {
    my $self = shift;
    return <<"EOS"
    genome processing-profile create clin-seq  --name 'November 2011 Clinical Sequencing'
    genome model define clin-seq  --processing-profile='November 2011 Clinical Sequencing'  --wgs-model=2882504846 --exome-model=2882505032 --tumor-rnaseq-model=2880794613
    # Automatically builds if/when the models have a complete underlying build
EOS
}

sub help_detail_for_create_profile {
    return <<EOS
The ClinSeq pipeline has the following processing-profile parameters.
  bam_readcount_version => version of bamreadcounts to use for read-counting steps,
  sireport_min_tumor_vaf => minimum tumor VAF cutoff to filter out variants in SnvIndelReport,
  sireport_max_normal_vaf   => maximum normal VAF cutoff to filter out variants in SnvIndelReport,
  sireport_min_coverage => minimum coverage cutoff to filter out variants in SnvIndelReport,
  sireport_min_mq_bq => semicolon delimited list of comma separated mapping quality and base quality scores, these
    are used as cutoffs for read-counting in SnvIndelReport.
EOS
}

sub _help_detail_for_model_define {
    return <<EOS

The ClinSeq pipeline takes four models, each of which is optional, and produces data sets potentially useful in a clinical or discovery project.

There are several primary goals:

1. Generate results that help to identify clinically actionable events (gain of function mutations, amplifications, etc.)

2. Converge results across data types for a single case (e.g. variant read counts from wgs + exome + rnaseq)

3. Automatically generate statistics, tables and figures to be used in manuscripts

4. Test novel analysis methods in a pipeline that runs faster and is more agile than say somatic-variation

EOS
}

sub _resolve_subject {
  my $self = shift;
  my @subjects = $self->_infer_candidate_subjects_from_input_models();
  if (@subjects > 1) {
    if ($self->force){
      @subjects = ($subjects[0]);
    }else{
      die $self->error_message(
          "Conflicting subjects on input models!:\n\t"
          . join("\n\t", map { $_->__display_name__ } @subjects)
      );
    }
  }elsif (@subjects == 0) {
     die $self->error_message("No subjects on input models?");
  }
  return $subjects[0];
}

sub _resolve_annotation {
  my $self = shift;
  my @annotations = $self->_infer_annotations_from_input_models();
  if (@annotations > 1) {
    if ($self->force){
      @annotations = ($annotations[0]);
    }else{
      die $self->error_message(
          "Conflicting annotations on input models!:\n\t"
          . join("\n\t", map { $_->__display_name__ } @annotations)
      );
    }
  }elsif (@annotations == 0) {
     die $self->error_message("No annotation builds on input models?");
  }
  return $annotations[0];
}

sub _resolve_reference {
  my $self = shift;
  my @references = $self->_infer_references_from_input_models();
  if (@references > 1) {
    if ($self->force){
      @references = ($references[0]);
    }else{
      die $self->error_message(
          "Conflicting reference sequence builds on input models!:\n\t"
          . join("\n\t", map { $_->__display_name__ } @references)
      );
    }
  }elsif (@references == 0) {
     die $self->error_message("No reference builds on input models?");
  }
  return $references[0];
}


#Implement specific error checking here, any error that is added to the @errors array will prevent the model from being commited to the database
#Could also implement a --force input above to allow over-riding of errors
sub __errors__ {
  my $self = shift;
  my @errors = $self->SUPER::__errors__;
  return @errors;
}

sub _initialize_build {
  my $self = shift;
  my $build = shift;

  # this is currently tracked as a metric on the build but it should really be an output/metric
  my $common_name = $self->expected_common_name;
  $build->common_name($common_name);

  return 1;
}

#MAP WORKFLOW INPUTS
sub map_workflow_inputs {
  my $self = shift;
  my $build = shift;

  my $model = $self;
  my $data_directory = $build->data_directory;

  # inputs
  my $wgs_build           = $build->wgs_build;
  my $exome_build         = $build->exome_build;
  my $tumor_rnaseq_build  = $build->tumor_rnaseq_build;
  my $normal_rnaseq_build = $build->normal_rnaseq_build;

  # set during build initialization
  my $common_name         = $build->common_name;
  unless ($common_name) {
      die "no common name set on build during initialization?";
  }

  my $cancer_annotation_db = $build->cancer_annotation_db;
  my $misc_annotation_db = $build->misc_annotation_db;
  my $cosmic_annotation_db = $build->cosmic_annotation_db;
  my ($mqs, $bqs) = $self->parse_qualities;

  # initial inputs used for various steps
  my @inputs = (
      model => $model,
      build => $build,
      wgs_build => $wgs_build,
      exome_build => $exome_build,
      tumor_rnaseq_build => $tumor_rnaseq_build,
      normal_rnaseq_build => $normal_rnaseq_build,
      working_dir => $data_directory,
      common_name => $common_name,
      verbose => 1,
      cancer_annotation_db => $cancer_annotation_db,
      misc_annotation_db => $misc_annotation_db,
      cosmic_annotation_db => $cosmic_annotation_db,
      bam_readcount_version => $self->bam_readcount_version,
      sireport_min_tumor_vaf => $self->sireport_min_tumor_vaf,
      sireport_max_normal_vaf => $self->sireport_max_normal_vaf,
      sireport_min_coverage => $self->sireport_min_coverage,
  );

  my $annotation_build = $self->_resolve_annotation;
  push @inputs, annotation_build => $annotation_build;

  #Identify test model
  if ($self->name =~ /^apipe\-test/){
      push @inputs, test => 1;
  }

  my $patient_dir = $data_directory . "/" . $common_name;
  my @dirs = ($patient_dir);

  #SummarizeBuilds
  my $input_summary_dir = $patient_dir . "/input";
  push @dirs, $input_summary_dir;
  push @inputs, summarize_builds_outdir => $input_summary_dir;
  push @inputs, summarize_builds_log_file => $input_summary_dir . "/QC_Report.tsv";

  if ($build->id) {
    push @inputs, summarize_builds_skip_lims_reports => 0;
  }else {
    #Watch out for -ve build IDs which will occur when the ClinSeq.t test is run.  In that case, do not run the LIMS reports
    push @inputs, summarize_builds_skip_lims_reports => 1;
  }

  #DumpIgvXml
  my $igv_session_dir = $patient_dir . '/igv';
  push @dirs, $igv_session_dir;
  push @inputs, igv_session_dir => $igv_session_dir;

  #GetVariantSources
  if ($wgs_build or $exome_build) {
    my $variant_sources_dir = $patient_dir . '/variant_source_callers';
    push @dirs, $variant_sources_dir;
    if ($wgs_build){
      my $wgs_variant_sources_dir = $variant_sources_dir . '/wgs';
      push @inputs, wgs_variant_sources_dir => $wgs_variant_sources_dir;
      push @dirs, $wgs_variant_sources_dir;
    }
    if ($exome_build){
      my $exome_variant_sources_dir = $variant_sources_dir . '/exome';
      push @inputs, exome_variant_sources_dir => $exome_variant_sources_dir;
      push @dirs, $exome_variant_sources_dir;
    }
  }

  #ImportSnvsIndels
  push @inputs, import_snvs_indels_filter_mt => 1;
  push @inputs, import_snvs_indels_outdir => $patient_dir;

  #CreateMutationDiagrams
  if ($wgs_build or $exome_build) {
    my $mutation_diagram_dir = $patient_dir . '/mutation_diagrams';
    push @dirs, $mutation_diagram_dir;
    push @inputs, (
        mutation_diagram_outdir => $mutation_diagram_dir,
        mutation_diagram_collapse_variants=>1,
        mutation_diagram_max_snvs_per_file=>1500,
        mutation_diagram_max_indels_per_file=>1500,
    );
  }

  #rnaseq analysis steps according to what rna-seq builds are defined as inputs
  if ($normal_rnaseq_build or $tumor_rnaseq_build){
    my $rnaseq_dir = $patient_dir . '/rnaseq';
    push @dirs, $rnaseq_dir;
    push @inputs, cufflinks_percent_cutoff => 1;

    #TophatJunctionsAbsolute and CufflinksExpressionAbsolute for 'normal' sample
    if ($normal_rnaseq_build){
      my $normal_rnaseq_dir = $rnaseq_dir . '/normal';
      push @dirs, $normal_rnaseq_dir;
      my $normal_tophat_junctions_absolute_dir = $normal_rnaseq_dir . '/tophat_junctions_absolute';
      push @dirs, $normal_tophat_junctions_absolute_dir;
      push @inputs, normal_tophat_junctions_absolute_dir => $normal_tophat_junctions_absolute_dir;
      my $normal_cufflinks_expression_absolute_dir = $normal_rnaseq_dir . '/cufflinks_expression_absolute';
      push @dirs, $normal_cufflinks_expression_absolute_dir;
      push @inputs, normal_cufflinks_expression_absolute_dir => $normal_cufflinks_expression_absolute_dir;
    }

    #TophatJunctionsAbsolute and CufflinksExpressionAbsolute for 'tumor' sample
    if ($tumor_rnaseq_build){
      my $tumor_rnaseq_dir = $rnaseq_dir . '/tumor';
      push @dirs, $tumor_rnaseq_dir;
      my $tumor_tophat_junctions_absolute_dir = $tumor_rnaseq_dir . '/tophat_junctions_absolute';
      push @dirs, $tumor_tophat_junctions_absolute_dir;
      push @inputs, tumor_tophat_junctions_absolute_dir => $tumor_tophat_junctions_absolute_dir;
      my $tumor_cufflinks_expression_absolute_dir = $tumor_rnaseq_dir . '/cufflinks_expression_absolute';
      push @dirs, $tumor_cufflinks_expression_absolute_dir;
      push @inputs, tumor_cufflinks_expression_absolute_dir => $tumor_cufflinks_expression_absolute_dir;
    }

    #CufflinksDifferentialExpression
    if ($normal_rnaseq_build and $tumor_rnaseq_build){
      my $cufflinks_differential_expression_dir = $rnaseq_dir . '/cufflinks_differential_expression';
      push @dirs, $cufflinks_differential_expression_dir;
      push @inputs, cufflinks_differential_expression_dir => $cufflinks_differential_expression_dir;
    }

    #Filtered and Intersected ChimeraScan fusion output for tumor RNAseq
    #The following will only work if the rna-seq build used a processing profile that uses chimerascan.
    if ($tumor_rnaseq_build){
      #Check for ChimeraScan fusion results
      if(-e $tumor_rnaseq_build->data_directory . '/fusions/filtered_chimeras.bedpe'){
        #copy over fusion files to this dir even if SV calls do not exist.
        my $tumor_filtered_fusion_dir = $patient_dir . '/rnaseq/tumor/fusions';
        push @dirs, $tumor_filtered_fusion_dir;
        if ($wgs_build){
          #Check for SV calls file
          if(-e $wgs_build->data_directory . '/effects/svs.hq.annotated'){
            my $ncbi_human_ensembl_build_id = $tumor_rnaseq_build->annotation_build->id;
            my $tumor_filtered_fusion_file =  $tumor_filtered_fusion_dir . '/filtered_chimeras.bedpe';
            my $wgs_sv_file = $build->wgs_build->data_directory . '/effects/svs.hq.annotated';
            my $tumor_filtered_intersected_fusion_file =  $tumor_filtered_fusion_dir . '/chimeras.filtered.intersected.bedpe';
            push @inputs, ncbi_human_ensembl_build_id => $ncbi_human_ensembl_build_id;
            push @inputs, wgs_sv_file => $wgs_sv_file;
            push @inputs, tumor_filtered_fusion_file => $tumor_filtered_fusion_file;
            push @inputs, tumor_filtered_intersected_fusion_file => $tumor_filtered_intersected_fusion_file;
          }
        }
      }
    }
  }

    #GenerateClonalityPlots
    if ($wgs_build){
      my $clonality_dir = $patient_dir . "/clonality/";
      push @dirs, $clonality_dir;
      push @inputs, clonality_dir => $clonality_dir;
    }

    #Make base-dir for CNV's
    if ($wgs_build or $exome_build or $self->has_microarray_build()){
      my $cnv_dir = $patient_dir . "/cnv/";
      push @dirs, $cnv_dir;
      push @inputs, cnv_dir => $cnv_dir;
    }

    #Make base-dir for WGS CNV's
    if ($wgs_build){
      my $wgs_cnv_dir = $patient_dir . "/cnv/wgs_cnv";
      push @dirs, $wgs_cnv_dir;
      push @inputs, wgs_cnv_dir => $wgs_cnv_dir;
    }

    #RunMicroArrayCnView
    if ($self->has_microarray_build()) {
      my $microarray_cnv_dir = $patient_dir . "/cnv/microarray_cnv/";
      push @dirs, $microarray_cnv_dir;
      push @inputs, microarray_cnv_dir => $microarray_cnv_dir;
    }

    #ExomeCNV
    if ($build->should_run_exome_cnv) {
      my $exome_cnv_dir = $patient_dir . "/cnv/exome_cnv/";
      push @dirs, $exome_cnv_dir;
      push @inputs, exome_cnv_dir => $exome_cnv_dir;
    }

    my $docm_variants_file = $self->_get_docm_variants_file($self->cancer_annotation_db);
    #Make DOCM report
    if (($exome_build or $wgs_build) and $docm_variants_file) {
      my $docm_report_dir = $patient_dir . "/docm_report/";
      push @dirs, $docm_report_dir;
      push @inputs, docm_report_dir => $docm_report_dir;
      push @inputs, docm_variants_file => $docm_variants_file;
      push @inputs, docmreport_min_bq => 0;
      push @inputs, docmreport_min_mq => 1;
    }

    #SummarizeSvs
    if ($wgs_build){
      my $sv_dir = $patient_dir . "/sv/";
      push @dirs, $sv_dir;
      push @inputs, sv_dir => $sv_dir;
    }

    #Make base-dir for WGS CNV's
    if ($wgs_build){
      my $wgs_cnv_summary_dir = $patient_dir . "/cnv/wgs_cnv/summary/";
      push @dirs, $wgs_cnv_summary_dir;
      push @inputs, wgs_cnv_summary_dir => $wgs_cnv_summary_dir;
    }

    #CreateMutationSpectrum
    my $iterator = List::MoreUtils::each_arrayref([1..@$mqs], $mqs, $bqs);
    while(my ($i, $mq, $bq) = $iterator->()) {
      if ($wgs_build) {
        push @inputs, 'wgs_mutation_spectrum_outdir' . $i =>
        $patient_dir . "/mutation-spectrum/b". $bq . "_q" . $mq .  "/";
        push @inputs, 'wgs_mutation_spectrum_datatype' . $i =>
        'wgs';
      }
      if ($exome_build) {
        push @inputs, 'exome_mutation_spectrum_outdir' . $i  =>
        $patient_dir . "/mutation-spectrum/b". $bq . "_q" . $mq .  "/";
        push @inputs, 'exome_mutation_spectrum_datatype' . $i =>
        'exome';
      }
    }
    if($wgs_build or $exome_build) {
      my $mutation_spectrum_dir = $patient_dir . "/mutation-spectrum";
      push @dirs, $mutation_spectrum_dir;
    }

    #AnnotateGenesByCategory
    push @inputs, 'gene_name_columns' => ['mapped_gene_name'];
    push @inputs, 'gene_name_regex' => 'mapped_gene_name';

    #MakeCircosPlot
    if ($wgs_build){
      my $circos_dir = $patient_dir . "/circos";
      push @dirs, $circos_dir;
      push @inputs, circos_outdir => $circos_dir;
    }

    #Converge SnvIndelReport and SciClone
    if ($exome_build || $wgs_build) {
      my $iterator = List::MoreUtils::each_arrayref([1..@$mqs], $mqs, $bqs);
      while(my ($i, $mq, $bq) = $iterator->()) {
        my $snv_indel_report_dir1 = File::Spec->join($patient_dir,
                                                     "/snv_indel_report/" . "b" . $bq . "_" . "q" . $mq);
        push @dirs, $snv_indel_report_dir1;
        push @inputs, "snv_indel_report_dir" . $i => $snv_indel_report_dir1;
        if($wgs_build or $build->should_run_exome_cnv) {
            my $sciclone_dir1 = File::Spec->join($patient_dir,
                                                 "/clonality/sciclone/" . "b" . $bq . "_" . "q" . $mq);
            push @dirs, $sciclone_dir1;
            push @inputs, "sciclone_dir" . $i  => $sciclone_dir1;
        }
        push @inputs, "sireport_min_bq" . $i => $bq;
        push @inputs, "sireport_min_mq" . $i => $mq;
      }
      my $target_gene_list = $cancer_annotation_db->data_directory . "/CancerGeneCensus/cancer_gene_census_ensgs.tsv";
      my $target_gene_list_name = "CancerGeneCensus";

      push @inputs, snv_indel_report_clean => 1;
      push @inputs, snv_indel_report_tmp_space => 1;
      push @inputs, snv_indel_report_target_gene_list => $target_gene_list;
      push @inputs, snv_indel_report_target_gene_list_name => $target_gene_list_name;
      push @inputs, snv_indel_report_tiers => 'tier1';
    }

    #IdentifyLoh
    if ($exome_build || $wgs_build) {
      my $loh_dir = $patient_dir . "/loh/";
      push @dirs, $loh_dir;
      push @inputs, loh_output_dir => $loh_dir;
    }

    # For now it works to create directories here because the data_directory has been allocated.
    #It is possible that this would not happen until later, which would mess up assigning inputs to many of the commands.
    for my $dir (@dirs) {
      Genome::Sys->create_directory($dir);
    }

    return @inputs;
}

sub _resolve_workflow_for_build {
  # This is called by Genome::Model::Build::start()
  # Returns a Workflow::Operation
  my $self = shift;
  my $build = shift;
  my $lsf_queue = shift;   # TODO: the workflow shouldn't need this yet
  my $lsf_project = shift;
  my ($mqs, $bqs) = $self->parse_qualities;
  if (!defined $lsf_queue || $lsf_queue eq '' || $lsf_queue eq 'inline') {
      $lsf_queue = Genome::Config::get('lsf_queue_build_worker_alt');
  }
  if (!defined $lsf_project || $lsf_project eq '') {
      $lsf_project = 'build' . $build->id;
  }

  # The wf system will call this method after we finish
  # but it consolidates logic to make a workflow which will
  # automatically take all inputs that method will assign.
  my %inputs = $self->map_workflow_inputs($build);
  my @input_properties = sort keys %inputs;

  # This must be updated for each new tool added which is "terminal" in the workflow!
  # (too bad it can't just be inferred from a dynamically expanding output connector)
  my @output_properties = qw(
      summarize_builds_result
      igv_session_result
  );

  if ($build->wgs_build) {
      push @output_properties, qw(
          summarize_wgs_tier1_snv_support_result
          summarize_svs_result
          summarize_cnvs_result
          clonality_result
          run_cn_view_result
          gene_category_cnv_amp_result
          gene_category_cnv_del_result
          gene_category_cnv_ampdel_result
          gene_category_wgs_snv_result
          gene_category_wgs_indel_result
          dgidb_cnv_amp_result
          dgidb_sv_fusion_result
          dgidb_wgs_snv_result
          dgidb_wgs_indel_result
          wgs_variant_sources_result
          circos_result
      );
  }

  if ($build->exome_build) {
      push @output_properties, qw(
          summarize_exome_tier1_snv_support_result
          gene_category_exome_snv_result
          gene_category_exome_indel_result
          dgidb_exome_snv_result
          dgidb_exome_indel_result
          exome_variant_sources_result
      );
  }

  if ($build->wgs_build and $build->exome_build) {
      push @output_properties, 'summarize_wgs_exome_tier1_snv_support_result';
      push @output_properties, 'gene_category_wgs_exome_indel_result';
      push @output_properties, 'gene_category_wgs_exome_snv_result';
      push @output_properties, 'dgidb_wgs_exome_indel_result';
      push @output_properties, 'dgidb_wgs_exome_snv_result';
  }

  if ($build->wgs_build or $build->exome_build) {
    push @output_properties, 'mutation_diagram_result';
    push @output_properties, 'import_snvs_indels_result';
    push @output_properties, 'loh_result';
    my $iterator = List::MoreUtils::each_arrayref([1..@$mqs], $mqs, $bqs);
    while(my ($i, $mq, $bq) = $iterator->()) {
      if($build->wgs_build) {
        push @output_properties,  'wgs_mutation_spectrum_result' . $i;
        push @output_properties, 'sciclone_result' . $i;
      }
      if($build->exome_build) {
        push @output_properties, 'exome_mutation_spectrum_result' . $i;
      }
      if($build->should_run_exome_cnv) {
        push @output_properties, 'sciclone_result' . $i;
      }
      push @output_properties, 'converge_snv_indel_report_result' . $i;
    }
  }

  if ($self->has_microarray_build()) {
      push @output_properties, 'microarray_cnv_result';
  }

  if ($build->should_run_exome_cnv) {
      push @output_properties, 'exome_cnv_result';
  }

  if (($build->exome_build or $build->wgs_build)
          and $self->_get_docm_variants_file($self->cancer_annotation_db)) {
      push @output_properties, 'docm_report_result';
  }

  if ($build->normal_rnaseq_build){
      push @output_properties, 'normal_tophat_junctions_absolute_result';
      push @output_properties, 'normal_cufflinks_expression_absolute_result';
  }

  if ($build->tumor_rnaseq_build){
      push @output_properties, 'tumor_tophat_junctions_absolute_result';
      push @output_properties, 'tumor_cufflinks_expression_absolute_result';
      push @output_properties, 'gene_category_cufflinks_result';
      push @output_properties, 'gene_category_tophat_result';
      push @output_properties, 'dgidb_cufflinks_result';
      push @output_properties, 'dgidb_tophat_result';
      if(-e $build->tumor_rnaseq_build->data_directory . '/fusions/filtered_chimeras.bedpe'){
          if ($build->wgs_build){
              if(-e $build->wgs_build->data_directory . '/effects/svs.hq.annotated'){
                  push @output_properties, 'intersect_tumor_fusion_sv_result';
              }
          }
      }
  }

  if ($build->normal_rnaseq_build and $build->tumor_rnaseq_build){
      push @output_properties, 'cufflinks_differential_expression_result';
      push @output_properties, 'gene_category_coding_de_up_result';
      push @output_properties, 'gene_category_coding_de_down_result';
      push @output_properties, 'gene_category_coding_de_result';
  }

  # Make the workflow and some convenience wrappers for adding steps and links

  my $workflow = Workflow::Model->create(
      name => $build->workflow_name,
      input_properties  => \@input_properties,
      output_properties => \@output_properties,
  );

  my $log_directory = $build->log_directory;
  $workflow->log_dir($log_directory);

  my $input_connector  = $workflow->get_input_connector;
  my $output_connector = $workflow->get_output_connector;

  my %steps_by_name;
  my $step = 0;
  my $add_step;
  my $add_link;

  $add_step = sub {
      my ($name, $cmd) = @_;

      if (substr($cmd,0,2) eq '::') {
          $cmd = 'Genome::Model::ClinSeq::Command' . $cmd;
      }

      unless ($cmd->can("execute")) {
          die "bad command $cmd!";
      }

      die "$name already used!" if $steps_by_name{$name};

      my $op = $workflow->add_operation(
          name => $name,
          operation_type => Workflow::OperationType::Command->create(
              command_class_name => $cmd,
          )
      );
      $op->operation_type->lsf_queue($lsf_queue);
      $op->operation_type->lsf_project($lsf_project);
      $steps_by_name{$name} = $op;

      # for inputs with these names, always add a link from the input connector
      # should this be default behavior for commands in a model's namespace?
      my $cmeta = $cmd->__meta__;
      for my $pname (qw/
          cancer_annotation_db
          misc_annotation_db
          cosmic_annotation_db
      /) {
          my $pmeta = $cmeta->property($pname);
          if ($pmeta) {
              $add_link->($input_connector, $pname, $op, $pname);
          }
      }

      return $op;
  };

  my $converge;

  $add_link = sub {
      my ($from_op, $from_p, $to_op, $to_p) = @_;
      $to_p = $from_p if not defined $to_p;

      if (ref($to_p) eq 'ARRAY') {
          Carp::confess("the 'to' property in a link cannot be a list!");
      }

      my $link;
      if (ref($from_p) eq 'ARRAY') {
          my $cname = "Combine: (@$from_p) for \"" . $to_op->name . "\" parameter \'$to_p\'";
          my $converge_op = $converge->($cname,$from_op,$from_p);
          $link = $workflow->add_link(
              left_operation => $converge_op,
              left_property => 'outputs',
              right_operation => $to_op,
              right_property => $to_p
          );
      }
      else {
          $link = $workflow->add_link(
              left_operation => $from_op,
              left_property => $from_p,
              right_operation => $to_op,
              right_property => $to_p
          );
      }
      $link or die "Failed to make link from $link from $from_p to $to_p!";
      return $link;
  };

  my $combo_cnt = 0;
  $converge = sub {
      my $cname = shift;
      my @params1 = @_;
      my @params2 = @_;

      # make a friendly name
      my $name = '';
      my $combo_cnt++;
      my $input_count = 0;
      while (@params1) {
          my $from_op = shift @params1;
          my $from_props = shift @params1;
          unless ($from_props) {
              die "expected \$op2,['p1','p2',...],\$op2,['p3','p4'],...";
          }
          unless (ref($from_props) eq 'ARRAY') {
              die "expected the second param (and every even param) in converge to be an arrayref of property names";
          }
          $input_count += scalar(@$from_props);
          if ($name) {
              $name .= " and ";
          }
          else {
              $name = "combine ($combo_cnt): ";
          }
          if ($from_op->name eq 'input_conector') {
              $name = "($combo_cnt)";
          }
          else {
              $name .= $from_op->name . "(";
          }
          for my $p (@$from_props) {
              $name .= $p;
              $name .= "," unless $p eq $from_props->[-1];
          }
          $name .= ")";
      }
      my $op = $workflow->add_operation(
          name => $cname,
          operation_type => Workflow::OperationType::Converge->create(
              input_properties => [ map { "i$_" } (1..$input_count) ],
              output_properties => ['outputs'],
          )
      );

      # create links
      my $input_n = 0;
      while (@params2) {
          my $from_op = shift @params2;
          my $from_props = shift @params2;
          for my $from_prop (@$from_props) {
              $input_n++;
              $add_link->(
                  $from_op,
                  $from_prop,
                  $op,
                  "i$input_n"
              );
          }
      }

      # return the op which will have a single "output"
      return $op;
  };

  #SummarizeBuilds - Summarize build inputs using SummarizeBuilds.pm
  my $msg = "Creating a summary of input builds using summarize-builds";
  my $summarize_builds_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeBuilds");
  $add_link->($input_connector, 'build', $summarize_builds_op, 'builds');
  $add_link->($input_connector, 'summarize_builds_outdir', $summarize_builds_op, 'outdir');
  $add_link->($input_connector, 'summarize_builds_skip_lims_reports', $summarize_builds_op, 'skip_lims_reports');
  $add_link->($input_connector, 'summarize_builds_log_file', $summarize_builds_op, 'log_file');
  $add_link->($summarize_builds_op, 'result', $output_connector, 'summarize_builds_result');

  #ImportSnvsIndels - Import SNVs and Indels
  my $import_snvs_indels_op;
  if ($build->wgs_build or $build->exome_build){
    my $msg = "Importing SNVs and Indels from somatic results, parsing, and merging exome/wgs if possible";
    $import_snvs_indels_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::ImportSnvsIndels");
    $add_link->($input_connector, 'import_snvs_indels_outdir', $import_snvs_indels_op, 'outdir');
    $add_link->($input_connector, 'import_snvs_indels_filter_mt', $import_snvs_indels_op, 'filter_mt');
    if ($build->wgs_build){
      $add_link->($input_connector, 'wgs_build', $import_snvs_indels_op, 'wgs_build');
    }
    if ($build->exome_build){
      $add_link->($input_connector, 'exome_build', $import_snvs_indels_op, 'exome_build');
    }
    $add_link->($import_snvs_indels_op, 'result', $output_connector, 'import_snvs_indels_result');
  }

  #GetVariantSources - Determine source variant caller for SNVs and InDels for wgs data
  my $wgs_variant_sources_op;
  if ($build->wgs_build) {
    my $msg = "Determining source variant callers of all tier1-3 SNVs and InDels for wgs data";
    $wgs_variant_sources_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::GetVariantSources");
    $add_link->($input_connector, 'wgs_build', $wgs_variant_sources_op, 'builds');
    $add_link->($input_connector, 'wgs_variant_sources_dir', $wgs_variant_sources_op, 'outdir');
    $add_link->($wgs_variant_sources_op, 'result', $output_connector, 'wgs_variant_sources_result');
  }
  #GetVariantSources - Determine source variant caller for SNVs and InDels for exome data
  my $exome_variant_sources_op;
  if ($build->exome_build) {
    my $msg = "Determining source variant callers of all tier1-3 SNVs and InDels for exome data";
    $exome_variant_sources_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::GetVariantSources");
    $add_link->($input_connector, 'exome_build', $exome_variant_sources_op, 'builds');
    $add_link->($input_connector, 'exome_variant_sources_dir', $exome_variant_sources_op, 'outdir');
    $add_link->($exome_variant_sources_op, 'result', $output_connector, 'exome_variant_sources_result');
  }

  #CreateMutationDiagrams - Create mutation diagrams (lolliplots) for all Tier1 SNVs/Indels and compare to COSMIC SNVs/Indels
  if ($build->wgs_build or $build->exome_build) {
    my $msg = "Creating mutation-diagram plots";
    my $mutation_diagram_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::CreateMutationDiagrams");
    if ($build->wgs_build and $build->exome_build) {
        $add_link->($input_connector, ['wgs_build','exome_build'], $mutation_diagram_op, 'builds');
    } elsif ($build->wgs_build) {
        $add_link->($input_connector, 'wgs_build', $mutation_diagram_op, 'builds');
    } elsif ($build->exome_build) {
        $add_link->($input_connector, 'exome_build', $mutation_diagram_op, 'builds');
    }
    $add_link->($mutation_diagram_op,'result',$output_connector,'mutation_diagram_result');

    for my $p (qw/outdir collapse_variants max_snvs_per_file max_indels_per_file/) {
        my $input_name = 'mutation_diagram_' . $p;
        $add_link->($input_connector,$input_name,$mutation_diagram_op,$p);
    }
  }

  #TophatJunctionsAbsolute - Run tophat junctions absolute analysis on normal
  my $normal_tophat_junctions_absolute_op;
  if ($build->normal_rnaseq_build){
    my $msg = "Performing tophat junction expression absolute analysis for normal sample";
    $normal_tophat_junctions_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::TophatJunctionsAbsolute');
    $add_link->($input_connector, 'normal_rnaseq_build', $normal_tophat_junctions_absolute_op, 'build');
    $add_link->($input_connector, 'normal_tophat_junctions_absolute_dir', $normal_tophat_junctions_absolute_op, 'outdir');
    $add_link->($normal_tophat_junctions_absolute_op, 'result', $output_connector, 'normal_tophat_junctions_absolute_result');
  }
  #TophatJunctionsAbsolute - Run tophat junctions absolute analysis on tumor
  my $tumor_tophat_junctions_absolute_op;
  if ($build->tumor_rnaseq_build){
    my $msg = "Performing tophat junction expression absolute analysis for tumor sample";
    $tumor_tophat_junctions_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::TophatJunctionsAbsolute');
    $add_link->($input_connector, 'tumor_rnaseq_build', $tumor_tophat_junctions_absolute_op, 'build');
    $add_link->($input_connector, 'tumor_tophat_junctions_absolute_dir', $tumor_tophat_junctions_absolute_op, 'outdir');
    $add_link->($tumor_tophat_junctions_absolute_op, 'result', $output_connector, 'tumor_tophat_junctions_absolute_result');
  }

  #CufflinksExpressionAbsolute - Run cufflinks expression absolute analysis on normal
  my $normal_cufflinks_expression_absolute_op;
  if ($build->normal_rnaseq_build){
    my $msg = "Performing cufflinks expression absolute analysis for normal sample";
    $normal_cufflinks_expression_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CufflinksExpressionAbsolute');
    $add_link->($input_connector, 'normal_rnaseq_build', $normal_cufflinks_expression_absolute_op, 'build');
    $add_link->($input_connector, 'normal_cufflinks_expression_absolute_dir', $normal_cufflinks_expression_absolute_op, 'outdir');
    $add_link->($input_connector, 'cufflinks_percent_cutoff', $normal_cufflinks_expression_absolute_op, 'percent_cutoff');
    $add_link->($normal_cufflinks_expression_absolute_op, 'result', $output_connector, 'normal_cufflinks_expression_absolute_result');
  }
  #CufflinksExpressionAbsolute - Run cufflinks expression absolute analysis on tumor
  my $tumor_cufflinks_expression_absolute_op;
  if ($build->tumor_rnaseq_build){
    my $msg = "Performing cufflinks expression absolute analysis for tumor sample";
    $tumor_cufflinks_expression_absolute_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CufflinksExpressionAbsolute');
    $add_link->($input_connector, 'tumor_rnaseq_build', $tumor_cufflinks_expression_absolute_op, 'build');
    $add_link->($input_connector, 'tumor_cufflinks_expression_absolute_dir', $tumor_cufflinks_expression_absolute_op, 'outdir');
    $add_link->($input_connector, 'cufflinks_percent_cutoff', $tumor_cufflinks_expression_absolute_op, 'percent_cutoff');
    $add_link->($tumor_cufflinks_expression_absolute_op, 'result', $output_connector, 'tumor_cufflinks_expression_absolute_result');
  }

  #CufflinksDifferentialExpression - Run cufflinks differential expression
  my $cufflinks_differential_expression_op;
  if ($build->normal_rnaseq_build and $build->tumor_rnaseq_build){
    my $msg = "Performing cufflinks differential expression analysis for case vs. control (e.g. tumor vs. normal)";
    $cufflinks_differential_expression_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CufflinksDifferentialExpression');
    $add_link->($input_connector, 'normal_rnaseq_build', $cufflinks_differential_expression_op, 'control_build');
    $add_link->($input_connector, 'tumor_rnaseq_build', $cufflinks_differential_expression_op, 'case_build');
    $add_link->($input_connector, 'cufflinks_differential_expression_dir', $cufflinks_differential_expression_op, 'outdir');
    $add_link->($cufflinks_differential_expression_op, 'result', $output_connector, 'cufflinks_differential_expression_result');
  }

  #Intersect filtered fusion calls with WGS SV calls.
  my $intersect_tumor_fusion_sv_op;
  if ($build->tumor_rnaseq_build){
    if(-e $build->tumor_rnaseq_build->data_directory . '/fusions/filtered_chimeras.bedpe'){
      #copy over fusion files
      $self->copy_fusion_files($build);
      if ($build->wgs_build){
        if(-e $build->wgs_build->data_directory . '/effects/svs.hq.annotated'){
          my $msg = "Intersecting filtered tumor ChimeraScan fusion calls with WGS SV calls.";
          $intersect_tumor_fusion_sv_op = $add_step->($msg, 'Genome::Model::Tools::ChimeraScan::IntersectSv');
          $add_link->($input_connector, 'ncbi_human_ensembl_build_id', $intersect_tumor_fusion_sv_op, 'annotation_build_id');
          $add_link->($input_connector, 'tumor_filtered_intersected_fusion_file', $intersect_tumor_fusion_sv_op, 'output_file');
          $add_link->($input_connector, 'wgs_sv_file', $intersect_tumor_fusion_sv_op, 'sv_output_file');
          $add_link->($input_connector, 'tumor_filtered_fusion_file', $intersect_tumor_fusion_sv_op, 'filtered_bedpe_file');
          $add_link->($intersect_tumor_fusion_sv_op, 'result', $output_connector, 'intersect_tumor_fusion_sv_result');
        }
      }
    }
  }

  #DumpIgvXml - Create IGV xml session files with increasing numbers of tracks and store in a single (WGS and Exome BAM files, RNA-seq BAM files, junctions.bed, SNV bed files, etc.)
  #genome model clin-seq dump-igv-xml --outdir=/gscuser/mgriffit/ --builds=119971814
  $msg = "Create IGV XML session files for varying levels of detail using the input builds";
  my $igv_session_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::DumpIgvXml");
  $add_link->($input_connector, 'build', $igv_session_op, 'builds');
  $add_link->($input_connector, 'igv_session_dir', $igv_session_op, 'outdir');
  $add_link->($igv_session_op, 'result', $output_connector, 'igv_session_result');

  #GenerateClonalityPlots - Run clonality analysis and produce clonality plots
  my $clonality_op;
  if ($build->wgs_build){
    my $msg = "Run clonality analysis and produce clonality plots";
    $clonality_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::GenerateClonalityPlots");
    $add_link->($input_connector, 'wgs_build', $clonality_op, 'somatic_var_build');
    $add_link->($input_connector, 'clonality_dir', $clonality_op, 'output_dir');
    $add_link->($input_connector, 'common_name', $clonality_op, 'common_name');
    $add_link->($input_connector, 'bam_readcount_version', $clonality_op, 'bam_readcount_version');
    $add_link->($clonality_op, 'result', $output_connector, 'clonality_result');
  }

  #RunCnView - Produce copy number results with run-cn-view.  Relies on clonality step already having been run
  my $run_cn_view_op;
  if ($build->wgs_build){
    my $msg = "Use gmt copy-number cn-view to produce copy number tables and images";
    $run_cn_view_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::RunCnView");
    $add_link->($input_connector, 'wgs_build', $run_cn_view_op, 'build');
    $add_link->($input_connector, 'wgs_cnv_dir', $run_cn_view_op, 'outdir');
    $add_link->($clonality_op, 'cnv_hmm_file', $run_cn_view_op);
    $add_link->($clonality_op, 'cnv_hq_file', $run_cn_view_op);
    $add_link->($run_cn_view_op, 'result', $output_connector, 'run_cn_view_result');
  }

  #RunMicroarrayCNV - produce cnv plots using microarray results
  my $microarray_cnv_op;
  if ($self->has_microarray_build()) {
    my $msg = "Call somatic copy number changes using microarray calls";
    $microarray_cnv_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::MicroarrayCnv");
    $add_link->($input_connector, 'microarray_cnv_dir', $microarray_cnv_op, 'outdir');
    $add_link->($input_connector, 'model', $microarray_cnv_op, 'clinseq_model');
    $add_link->($input_connector, 'annotation_build', $microarray_cnv_op, 'annotation_build_id');
    $add_link->($microarray_cnv_op, 'result', $output_connector, 'microarray_cnv_result');
  }

  #RunExomeCNV - produce cnv plots using WEx results
  my $exome_cnv_op;
  if ($build->should_run_exome_cnv) {
    my $msg = "Call somatic copy number changes using exome data";
    $exome_cnv_op = $add_step->($msg, "Genome::Model::Tools::CopyNumber::Cnmops");
    $add_link->($input_connector, 'exome_cnv_dir', $exome_cnv_op, 'outdir');
    $add_link->($input_connector, 'model', $exome_cnv_op, 'clinseq_model');
    $add_link->($input_connector, 'annotation_build', $exome_cnv_op, 'annotation_build_id');
    $add_link->($exome_cnv_op, 'result', $output_connector, 'exome_cnv_result');
  }

  #RunDOCMReport
  my $docm_report_op;
  if (($build->exome_build or $build->wgs_build)
        and $self->_get_docm_variants_file($self->cancer_annotation_db)) {
    my $msg = "Produce a report using DOCM";
    $docm_report_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::Converge::DocmReport");
    $add_link->($input_connector, 'docm_report_dir', $docm_report_op, 'outdir');
    $add_link->($input_connector, 'build', $docm_report_op, 'builds');
    $add_link->($input_connector, 'docm_variants_file', $docm_report_op, 'docm_variants_file');
    $add_link->($input_connector, 'bam_readcount_version', $docm_report_op, 'bam_readcount_version');
    $add_link->($input_connector, 'docmreport_min_bq', $docm_report_op, 'bq');
    $add_link->($input_connector, 'docmreport_min_mq', $docm_report_op, 'mq');
    $add_link->($docm_report_op, 'result', $output_connector, 'docm_report_result');
}

  #SummarizeCnvs - Generate a summary of CNV results, copy cnvs.hq, cnvs.png, single-bam copy number plot PDF, etc. to the cnv directory
  #This step relies on the generate-clonality-plots step already having been run
  #It also relies on run-cn-view step having been run already
  my $summarize_cnvs_op;
  if ($build->wgs_build){
    my $msg = "Summarize CNV results from WGS somatic variation";
    $summarize_cnvs_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeCnvs");
    $add_link->($input_connector, 'wgs_cnv_summary_dir', $summarize_cnvs_op, 'outdir');
    $add_link->($input_connector, 'wgs_build', $summarize_cnvs_op, 'build');
    $add_link->($clonality_op, 'cnv_hmm_file', $summarize_cnvs_op);
    $add_link->($clonality_op, 'cnv_hq_file', $summarize_cnvs_op);
    $add_link->($run_cn_view_op, 'gene_amp_file', $summarize_cnvs_op);
    $add_link->($run_cn_view_op, 'gene_del_file', $summarize_cnvs_op);
    $add_link->($summarize_cnvs_op, 'result', $output_connector, 'summarize_cnvs_result');
  }

  #SummarizeSvs - Generate a summary of SV results from the WGS SV results
  my $summarize_svs_op;
  if ($build->wgs_build){
    my $msg = "Summarize SV results from WGS somatic variation";
    $summarize_svs_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeSvs");
    $add_link->($input_connector, 'wgs_build', $summarize_svs_op, 'builds');
    $add_link->($input_connector, 'sv_dir', $summarize_svs_op, 'outdir');
    $add_link->($summarize_svs_op, 'result', $output_connector, 'summarize_svs_result');
  }

  #Add gene category annotations to some output files from steps above. (e.g. determine which SNV affected genes are kinases, ion channels, etc.)
  #AnnotateGenesByCategory - gene_category_exome_snv_result
  if ($build->exome_build){
    my $msg = "Add gene category annotations to SNVs identified by exome";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($import_snvs_indels_op, 'exome_snv_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_exome_snv_result');
  }
  #AnnotateGenesByCategory - gene_category_wgs_snv_result
  if ($build->wgs_build){
    my $msg = "Add gene category annotations to SNVs identified by wgs";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($import_snvs_indels_op, 'wgs_snv_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_snv_result');
  }
  #AnnotateGenesByCategory - gene_category_wgs_exome_snv_result
  if ($build->wgs_build and $build->exome_build){
    my $msg = "Add gene category annotations to SNVs identified by wgs OR exome";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($import_snvs_indels_op, 'wgs_exome_snv_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_exome_snv_result');
  }
  #AnnotateGenesByCategory - gene_category_exome_indel_result
  if ($build->exome_build){
    my $msg = "Add gene category annotations to InDels identified by exome";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($import_snvs_indels_op, 'exome_indel_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_exome_indel_result');
  }
  #AnnotateGenesByCategory - gene_category_wgs_indel_result
  if ($build->wgs_build){
    my $msg = "Add gene category annotations to InDels identified by wgs";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($import_snvs_indels_op, 'wgs_indel_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_indel_result');
  }
  #AnnotateGenesByCategory - gene_category_wgs_exome_indel_result
  if ($build->wgs_build and $build->exome_build){
    my $msg = "Add gene category annotations to InDels identified by wgs OR exome";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($import_snvs_indels_op, 'wgs_exome_indel_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_wgs_exome_indel_result');
  }
  #AnnotateGenesByCategory - gene_category_cnv_amp_result
  if ($build->wgs_build){
    my $msg = "Add gene category annotations to CNV amp genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($run_cn_view_op, 'gene_amp_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cnv_amp_result');
  }
  #AnnotateGenesByCategory - gene_category_cnv_del_result
  if ($build->wgs_build){
    my $msg = "Add gene category annotations to CNV del genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($run_cn_view_op, 'gene_del_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cnv_del_result');
  }
  #AnnotateGenesByCategory - gene_category_cnv_ampdel_result
  if ($build->wgs_build){
    my $msg = "Add gene category annotations to CNV amp OR del genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($run_cn_view_op, 'gene_ampdel_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cnv_ampdel_result');
  }
  #AnnotateGenesByCategory - gene_category_cufflinks_result
  if ($build->tumor_rnaseq_build){
    my $msg = "Add gene category annotations to cufflinks top1 percent genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($tumor_cufflinks_expression_absolute_op, 'tumor_fpkm_topnpercent_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_cufflinks_result');
  }
  #AnnotateGenesByCategory - gene_category_tophat_result
  if ($build->tumor_rnaseq_build){
    my $msg = "Add gene category annotations to tophat junctions top1 percent genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($tumor_tophat_junctions_absolute_op, 'junction_topnpercent_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_tophat_result');
  }
  #AnnotateGenesByCategory - gene_category_coding_de_up_result
  if ($build->normal_rnaseq_build && $build->tumor_rnaseq_build){
    my $msg = "Add gene category annotations to up-regulated coding genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($cufflinks_differential_expression_op, 'coding_hq_up_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_coding_de_up_result');
  }
  #AnnotateGenesByCategory - gene_category_coding_de_down_result
  if ($build->normal_rnaseq_build && $build->tumor_rnaseq_build){
    my $msg = "Add gene category annotations to down-regulated coding genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($cufflinks_differential_expression_op, 'coding_hq_down_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_coding_de_down_result');
  }
  #AnnotateGenesByCategory - gene_category_coding_de_result
  if ($build->normal_rnaseq_build && $build->tumor_rnaseq_build){
    my $msg = "Add gene category annotations to DE coding genes";
    my $annotate_genes_by_category_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByCategory");
    $add_link->($cufflinks_differential_expression_op, 'coding_hq_de_file', $annotate_genes_by_category_op, 'infile');
    $add_link->($input_connector, 'gene_name_columns', $annotate_genes_by_category_op, 'gene_name_columns');
    $add_link->($annotate_genes_by_category_op, 'result', $output_connector, 'gene_category_coding_de_result');
  }

  #DGIDB gene annotation
  if ($build->exome_build){
    $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'exome_snv_file', $input_connector, $output_connector, 'dgidb_exome_snv_result');
    $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'exome_indel_file', $input_connector, $output_connector, 'dgidb_exome_indel_result');
  }

  if ($build->wgs_build){
    $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_snv_file', $input_connector, $output_connector, 'dgidb_wgs_snv_result');
    $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_indel_file', $input_connector, $output_connector, 'dgidb_wgs_indel_result');
    $self->add_dgidb_op_to_flow($add_step, $add_link, $run_cn_view_op, 'gene_amp_file', $input_connector, $output_connector, 'dgidb_cnv_amp_result');
    $self->add_dgidb_op_to_flow($add_step, $add_link, $summarize_svs_op, 'fusion_output_file', $input_connector, $output_connector, 'dgidb_sv_fusion_result');
  }

  if ($build->wgs_build and $build->exome_build) {
    $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_exome_snv_file', $input_connector, $output_connector, 'dgidb_wgs_exome_snv_result');
    $self->add_dgidb_op_to_flow($add_step, $add_link, $import_snvs_indels_op, 'wgs_exome_indel_file', $input_connector, $output_connector, 'dgidb_wgs_exome_indel_result');
  }

  if ($build->tumor_rnaseq_build){
    $self->add_dgidb_op_to_flow($add_step, $add_link, $tumor_cufflinks_expression_absolute_op, 'tumor_fpkm_topnpercent_file', $input_connector, $output_connector, 'dgidb_cufflinks_result');
    $self->add_dgidb_op_to_flow($add_step, $add_link, $tumor_tophat_junctions_absolute_op, 'junction_topnpercent_file', $input_connector, $output_connector, 'dgidb_tophat_result');
  }

  #SummarizeTier1SnvSupport - For each of the following: WGS SNVs, Exome SNVs, and WGS+Exome SNVs, do the following:
  #Get BAM readcounts for WGS (tumor/normal), Exome (tumor/normal), RNAseq (tumor), RNAseq (normal) - as available of course
  #TODO: Break this down to do direct calls to GetBamReadCounts instead of wrapping it.
  my $summarize_tier1_snv_support_op;
  for my $run (qw/wgs exome wgs_exome/) {
    if ($run eq 'wgs' and not $build->wgs_build) {
      next;
    }
    if ($run eq 'exome' and not $build->exome_build) {
      next;
    }
    if ($run eq 'wgs_exome' and not ($build->wgs_build and $build->exome_build)) {
      next;
    }
    my $txt_name = $run;
    $txt_name =~ s/_/ plus /g;
    $txt_name =~ s/wgs/WGS/;
    $txt_name =~ s/exome/Exome/;
    $msg = "$txt_name Summarize Tier 1 SNV Support (BAM read counts)";
    $summarize_tier1_snv_support_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::SummarizeTier1SnvSupport");
    $add_link->($import_snvs_indels_op, $run . "_snv_file", $summarize_tier1_snv_support_op, $run . "_positions_file");
    $add_link->($input_connector, 'wgs_build', $summarize_tier1_snv_support_op);
    $add_link->($input_connector, 'exome_build', $summarize_tier1_snv_support_op);
    $add_link->($input_connector, 'tumor_rnaseq_build', $summarize_tier1_snv_support_op);
    $add_link->($input_connector, 'normal_rnaseq_build', $summarize_tier1_snv_support_op);
    if ($build->tumor_rnaseq_build){
      $add_link->($tumor_cufflinks_expression_absolute_op, 'tumor_fpkm_file', $summarize_tier1_snv_support_op);
    }
    $add_link->($input_connector, 'verbose', $summarize_tier1_snv_support_op);
    $add_link->($summarize_tier1_snv_support_op, 'result', $output_connector, "summarize_${run}_tier1_snv_support_result");
  }

  #MakeCircosPlot - Creates a Circos plot to summarize the data using MakeCircosPlot.pm
  #Currently WGS data is a minimum requirement for Circos plot generation.
  my $make_circos_plot_op;
  if ($build->wgs_build){
    $msg = "Create a Circos plot using MakeCircosPlot";
    $make_circos_plot_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::MakeCircosPlot");
    $add_link->($input_connector, 'build', $make_circos_plot_op);
    $add_link->($input_connector, 'circos_outdir', $make_circos_plot_op, 'output_directory');
    $add_link->($summarize_svs_op, 'fusion_output_file', $make_circos_plot_op, 'candidate_fusion_infile');
    $add_link->($clonality_op, 'cnv_hmm_file', $make_circos_plot_op);
    if($build->normal_rnaseq_build || $build->tumor_rnaseq_build){
      if($build->normal_rnaseq_build){
          $add_link->($cufflinks_differential_expression_op, 'coding_hq_de_file', $make_circos_plot_op);
      }else{
          $add_link->($tumor_cufflinks_expression_absolute_op, 'tumor_fpkm_topnpercent_file', $make_circos_plot_op);
      }
    }
    $add_link->($import_snvs_indels_op, 'result', $make_circos_plot_op, 'import_snvs_indels_result');
    $add_link->($run_cn_view_op, 'gene_ampdel_file', $make_circos_plot_op);
    $add_link->($make_circos_plot_op, 'result', $output_connector, 'circos_result');
  }

  #Converge SnvIndelReport
  my $converge_snv_indel_report_op1;
  my @converge_snv_indel_report_ops;
  if ($build->wgs_build || $build->exome_build) {
    #Create a report for each $bq $mq combo.
    my $iterator = List::MoreUtils::each_arrayref([1..@$mqs], $mqs, $bqs);
    while(my ($i, $mq, $bq) = $iterator->()) {
        $msg = "Generate SnvIndel Report" . $i . ".";
        $converge_snv_indel_report_op1 = $add_step->($msg, "Genome::Model::ClinSeq::Command::Converge::SnvIndelReport");
        $add_link->($input_connector, 'build', $converge_snv_indel_report_op1, 'builds');
        $add_link->($input_connector, 'bam_readcount_version', $converge_snv_indel_report_op1, 'bam_readcount_version');
        $add_link->($input_connector, 'snv_indel_report_dir' . $i, $converge_snv_indel_report_op1, 'outdir');
        $add_link->($input_connector, 'snv_indel_report_clean', $converge_snv_indel_report_op1, 'clean');
        $add_link->($input_connector, 'snv_indel_report_tmp_space', $converge_snv_indel_report_op1, 'tmp_space');
        $add_link->($input_connector, 'annotation_build', $converge_snv_indel_report_op1, 'annotation_build');
        $add_link->($input_connector, 'snv_indel_report_target_gene_list', $converge_snv_indel_report_op1, 'target_gene_list');
        $add_link->($input_connector, 'snv_indel_report_target_gene_list_name', $converge_snv_indel_report_op1, 'target_gene_list_name');
        $add_link->($input_connector, 'sireport_min_tumor_vaf', $converge_snv_indel_report_op1, 'min_tumor_vaf');
        $add_link->($input_connector, 'sireport_max_normal_vaf', $converge_snv_indel_report_op1, 'max_normal_vaf');
        $add_link->($input_connector, 'sireport_min_coverage', $converge_snv_indel_report_op1, 'min_coverage');
        $add_link->($input_connector, 'sireport_min_bq' . $i, $converge_snv_indel_report_op1, 'bq');
        $add_link->($input_connector, 'sireport_min_mq' . $i, $converge_snv_indel_report_op1, 'mq');
        if ($build->wgs_build){
            $add_link->($wgs_variant_sources_op, 'snv_variant_sources_file', $converge_snv_indel_report_op1, '_wgs_snv_variant_sources_file');
            $add_link->($wgs_variant_sources_op, 'indel_variant_sources_file', $converge_snv_indel_report_op1, '_wgs_indel_variant_sources_file');
        }
        if ($build->exome_build){
            $add_link->($exome_variant_sources_op, 'snv_variant_sources_file', $converge_snv_indel_report_op1, '_exome_snv_variant_sources_file');
            $add_link->($exome_variant_sources_op, 'indel_variant_sources_file', $converge_snv_indel_report_op1, '_exome_indel_variant_sources_file');
        }
        #If this is a build of a test model, perform a faster analysis (e.g. apipe-test-clinseq-wer)
        if ($self->name =~ /^apipe\-test/){
          $add_link->($input_connector, 'snv_indel_report_tiers', $converge_snv_indel_report_op1, 'tiers');
        }
        $add_link->($converge_snv_indel_report_op1, 'result', $output_connector, 'converge_snv_indel_report_result' . $i);
        push @converge_snv_indel_report_ops, $converge_snv_indel_report_op1;
    }
  }

  #CreateMutationDiagrams - Create mutation spectrum results for wgs data
  if ($build->wgs_build) {
    my $iterator = List::MoreUtils::each_arrayref([1..@$mqs], $mqs, $bqs);
    while(my ($i, $mq, $bq) = $iterator->()) {
      my $msg = "Creating mutation spectrum results for wgs snvs using create-mutation-spectrum" . $i;
      my $create_mutation_spectrum_wgs_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CreateMutationSpectrum');
      $add_link->($input_connector, 'build', $create_mutation_spectrum_wgs_op, 'clinseq_build');
      $add_link->($input_connector, 'wgs_build', $create_mutation_spectrum_wgs_op, 'somvar_build');
      $add_link->($input_connector, 'wgs_mutation_spectrum_outdir' . $i, $create_mutation_spectrum_wgs_op, 'outdir');
      $add_link->($input_connector, 'wgs_mutation_spectrum_datatype' . $i, $create_mutation_spectrum_wgs_op, 'datatype');
      $add_link->($input_connector, 'sireport_min_bq' . $i, $create_mutation_spectrum_wgs_op, 'min_base_quality');
      $add_link->($input_connector, 'sireport_min_mq' . $i, $create_mutation_spectrum_wgs_op, 'min_quality_score');
      $add_link->($converge_snv_indel_report_ops[$i-1], 'result', $create_mutation_spectrum_wgs_op, 'converge_snv_indel_report_result');
      $add_link->($create_mutation_spectrum_wgs_op, 'result', $output_connector, 'wgs_mutation_spectrum_result' . $i);
    }
  }

  #CreateMutationDiagrams - Create mutation spectrum results for exome data
  if ($build->exome_build) {
    my $iterator = List::MoreUtils::each_arrayref([1..@$mqs], $mqs, $bqs);
    while(my ($i, $mq, $bq) = $iterator->()) {
      my $msg = "Creating mutation spectrum results for exome snvs using create-mutation-spectrum " . $i;
      my $create_mutation_spectrum_exome_op = $add_step->($msg, 'Genome::Model::ClinSeq::Command::CreateMutationSpectrum');
      $add_link->($input_connector, 'build', $create_mutation_spectrum_exome_op, 'clinseq_build');
      $add_link->($input_connector, 'exome_build', $create_mutation_spectrum_exome_op, 'somvar_build');
      $add_link->($input_connector, 'exome_mutation_spectrum_outdir' . $i, $create_mutation_spectrum_exome_op, 'outdir' );
      $add_link->($input_connector, 'exome_mutation_spectrum_datatype' . $i, $create_mutation_spectrum_exome_op, 'datatype');
      $add_link->($input_connector, 'sireport_min_bq' . $i, $create_mutation_spectrum_exome_op, 'min_base_quality');
      $add_link->($input_connector, 'sireport_min_mq' . $i, $create_mutation_spectrum_exome_op, 'min_quality_score');
      $add_link->($converge_snv_indel_report_ops[$i-1], 'result', $create_mutation_spectrum_exome_op, 'converge_snv_indel_report_result');
      $add_link->($create_mutation_spectrum_exome_op, 'result', $output_connector, 'exome_mutation_spectrum_result' . $i);
    }
  }

  #GenerateSciClonePlots - Run clonality analysis and produce clonality plots
  if ($build->wgs_build or $build->should_run_exome_cnv){
    my $iterator = List::MoreUtils::each_arrayref([1..@$mqs], $mqs, $bqs);
    while(my ($i, $mq, $bq) = $iterator->()) {
      my $msg = "Run clonality analysis and produce clonality plots using SciClone " . $i;
      my $sciclone_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::GenerateSciclonePlots");
      $add_link->($input_connector, 'sciclone_dir' . $i, $sciclone_op, 'outdir');
      $add_link->($input_connector, 'build', $sciclone_op, 'clinseq_build');
      $add_link->($input_connector, 'sireport_min_coverage', $sciclone_op, 'min_coverage');
      $add_link->($input_connector, 'sireport_min_mq' . $i, $sciclone_op, 'min_mq');
      $add_link->($input_connector, 'sireport_min_bq' . $i, $sciclone_op, 'min_bq');
      if ($build->wgs_build) {
        $add_link->($run_cn_view_op, 'result', $sciclone_op, 'wgs_cnv_result');
      }
      if ($build->should_run_exome_cnv) {
        $add_link->($exome_cnv_op, 'result', $sciclone_op, 'exome_cnv_result');
      }
      if ($self->has_microarray_build()) {
        $add_link->($microarray_cnv_op, 'result', $sciclone_op, 'microarray_cnv_result');
      }
      $add_link->($converge_snv_indel_report_ops[$i-1], 'result', $sciclone_op, 'converge_snv_indel_report_result');
      $add_link->($sciclone_op, 'result', $output_connector, 'sciclone_result' . $i);
    }
  }

  #IdentifyLoh - Run identify-loh tool for exome or WGS data
  #genome model clin-seq identify-loh --clinseq-build=fafd219665d54462893fbacfe6639f70 --outdir=/Documents/GTB11/ --bamrc-version=0.7
  if ($build->wgs_build or $build->exome_build){
    $msg = "Identify regions of LOH and create plots";
    my $identify_loh_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::IdentifyLoh");
    $add_link->($input_connector, 'build', $identify_loh_op, 'clinseq_build');
    $add_link->($input_connector, 'loh_output_dir', $identify_loh_op, 'outdir');
    $add_link->($input_connector, 'bam_readcount_version', $identify_loh_op, 'bamrc_version');
    if ($self->name =~ /^apipe\-test/){
      $add_link->($input_connector, 'test', $identify_loh_op, 'test');
    }
    $add_link->($identify_loh_op, 'result', $output_connector, 'loh_result');
  }

  # REMINDER:
  # For new steps be sure to add their result to the output connector if they do not feed into another step.
  # When you do that, expand the list of output properties above.

  my @errors = $workflow->validate();
  if (@errors) {
      for my $error (@errors) {
          $self->error_message($error);
      }
      die "Invalid workflow!";
  }

  return $workflow;
}

sub add_dgidb_op_to_flow {
  my ($self, $add_step, $add_link, $op, $op_prop, $input_connector, $output_connector, $out_prop) = @_;
  my $msg = "Add dgidb gene annotations to $op_prop";
  my $annotate_genes_by_dgidb_op = $add_step->($msg, "Genome::Model::ClinSeq::Command::AnnotateGenesByDgidb");
  $add_link->($op, $op_prop, $annotate_genes_by_dgidb_op, 'input_file');
  $add_link->($input_connector, 'gene_name_regex', $annotate_genes_by_dgidb_op, 'gene_name_regex');
  $add_link->($annotate_genes_by_dgidb_op, 'result', $output_connector, $out_prop);
  return 1;
}


sub _infer_candidate_subjects_from_input_models {
  my $self = shift;
  
  my %subjects;
  foreach my $input_model ($self->_input_models){
    next unless $input_model;
    my $patient;
    if ($input_model->subject->isa("Genome::Individual")){
      $patient = $input_model->subject;
    } else {
      $patient = $input_model->subject->individual;
    }
    $subjects{ $patient->id } = $patient;

    if ($input_model->can("tumor_model")) {
      $subjects{ $input_model->tumor_model->subject->individual->id } = $input_model->tumor_model->subject->individual;
    }
    if ($input_model->can("normal_model")) {
      $subjects{ $input_model->normal_model->subject->individual->id } = $input_model->normal_model->subject->individual;
    }
  }
  my @subjects = sort { $a->id cmp $b->id } values %subjects;
  return @subjects;
}

sub _infer_annotations_from_input_models {
  my $self = shift;
  
  my %annotations;
  foreach my $input_model ($self->_input_models){
    next unless $input_model;
    if ($input_model->can("annotation_build")){
      my $annotation_build = $input_model->annotation_build;
      $annotations{ $annotation_build->id } = $annotation_build;
    }
    if ($input_model->can("annotation_reference_build")){
      my $annotation_build = $input_model->annotation_reference_build;
      $annotations{ $annotation_build->id } = $annotation_build;
    }
  }
  my @annotations = sort { $a->id cmp $b->id } values %annotations;
  return @annotations;
}

sub _infer_references_from_input_models {
  my $self = shift;
  
  my %references;
  foreach my $input_model ($self->_input_models){
    next unless $input_model;
    if ($input_model->can("reference_sequence_build")){
      my $reference_build = $input_model->reference_sequence_build;
      $references{ $reference_build->id } = $reference_build;
    }
  }
  my @references = sort { $a->id cmp $b->id } values %references;
  return @references;
}

sub _infer_input_models_for_somatic_model {
    my $self = shift;
    my $input_models = shift;
    my $somatic_model = shift;

    if ($somatic_model->class eq 'Genome::Model::SomaticValidation') { return; }

    push (@$input_models, $somatic_model->normal_model);
    push (@$input_models, $somatic_model->tumor_model);
}

sub _input_models {
    my $self = shift;
    
    my @input_models = ( $self->wgs_model, $self->exome_model, $self->tumor_rnaseq_model, $self->normal_rnaseq_model );

    $self->_infer_input_models_for_somatic_model(\@input_models,$self->wgs_model) if $self->wgs_model;
    $self->_infer_input_models_for_somatic_model(\@input_models,$self->exome_model) if $self->exome_model;

    return @input_models;
}

sub _resolve_resource_requirements_for_build {
  #Set LSF parameters for the ClinSeq job
  my $lsf_resource_string = "-R 'rusage[tmp=10000:mem=4000]' -M 4000000";
  return($lsf_resource_string);
}

# This is implemented here until refactoring is done on the model/build API
# It ensures that when the CI server compares two clinseq builds it ignores the correct files.
# Keep it in sync with the diff conditions in ClinSeq.t.
sub files_ignored_by_build_diff {
    return qw(
        build.xml
        reports/Build_Initialized/report.xml
        reports/Build_Succeeded/report.xml
        logs/.*
        .*.R$
        .*.Robject$
        .*.pdf$
        .*.jpg$
        .*.jpeg$
        .*.png$
        .*.svg$
        .*.xls$
        .*.bam$
        .*/summary/rc_summary.stderr$
        .*._COSMIC.svg$
        .*.clustered.data.tsv$
        .*QC_Report.tsv$
        .*.DumpIgvXml.log.txt$
        .*/mutation_diagrams/cosmic.mutation-diagram.stderr$
        .*/mutation_diagrams/somatic.mutation-diagram.stderr$
        .*.sequence-context.stderr$
        .*.annotate.stderr$
        .*/mutation-spectrum/exome/summarize_mutation_spectrum/summarize-mutation-spectrum.stderr$
        .*/mutation-spectrum/wgs/summarize_mutation_spectrum/summarize-mutation-spectrum.stderr$
        .*LIMS_Sample_Sequence_QC_library.tsv$
        .*.R.stderr$
        .*/snv_indel_report/subjects_legend.txt$
        .*sciclone.*.clusters.txt$
        .*mutation_spectrum.input.tsv$
        .*/kinase_only.tsv$
        .*/all_interactions.tsv$
        .*/expert_antineoplastic.tsv$
    );
};

# Below are some methods to retrieve files from a build
# Ideally they would be tied to creation of these file paths, but they currently
# throw an exception if the files don't exist. Consider another approach...
sub patient_dir {
  my ($self, $build) = @_;
  unless($build->common_name) {
    die $self->error_message("Common name is not defined.");
  }
  my $patient_dir = $build->data_directory . "/" . $build->common_name;
  unless(-d $patient_dir) {
    die $self->error_message("ClinSeq patient directory not found. Expected: $patient_dir");
  }
  return $patient_dir;
}

sub snv_variant_source_file {
  my ($self, $build, $data_type,) = @_;

  my $patient_dir = $self->patient_dir($build);
  my $source = $patient_dir . "/variant_source_callers";
  my $exception;
  if(-d $source) {
     my $dir = $source . "/$data_type";
     if(-d $dir) {
         my $file = $dir . "/snv_sources.tsv";
         if(-e $file) {
             return $file;
         }
         else {
             $exception = $self->error_message("Expected $file inside $dir and it did not exist.");
         }
     }
     else {
         $exception = $self->error_message("$data_type sub-directory not found in $source.");
     }
  }
  else {
      $exception = $self->error_message("$source directory not found");
  }
  die $exception;
}

sub copy_fusion_files {
  my ($self, $build) = @_;
  my $rnaseq_build_dir = $build->tumor_rnaseq_build->data_directory;
  my $tumor_unfiltered_fusion_file =  $rnaseq_build_dir . '/fusions/Genome_Model_RnaSeq_DetectFusionsResult_Chimerascan_VariableReadLength_Result/chimeras.bedpe';
  my $tumor_filtered_fusion_file =  $rnaseq_build_dir . '/fusions/filtered_chimeras.bedpe';
  my $tumor_filtered_annotated_fusion_file =  $rnaseq_build_dir . '/fusions/filtered_chimeras.catanno.bedpe';
  my $clinseq_tumor_unfiltered_fusion_file = $self->patient_dir($build) . '/rnaseq/tumor/fusions/chimeras.bedpe';
  my $clinseq_tumor_filtered_fusion_file = $self->patient_dir($build) . '/rnaseq/tumor/fusions/filtered_chimeras.bedpe';
  my $clinseq_tumor_filtered_annotated_fusion_file = $self->patient_dir($build) . '/rnaseq/tumor/fusions/filtered_chimeras.catanno.bedpe';
  if(-e $tumor_unfiltered_fusion_file) {
      unless(Genome::Sys->copy_file($tumor_unfiltered_fusion_file, $clinseq_tumor_unfiltered_fusion_file)) {
         die "unable to copy $tumor_unfiltered_fusion_file";
      }
  }
  if(-e $tumor_filtered_fusion_file) {
      unless(Genome::Sys->copy_file($tumor_filtered_fusion_file, $clinseq_tumor_filtered_fusion_file)) {
         die "unable to copy $tumor_filtered_fusion_file";
      }
  }
  if(-e $tumor_filtered_annotated_fusion_file) {
      unless(Genome::Sys->copy_file($tumor_filtered_annotated_fusion_file, $clinseq_tumor_filtered_annotated_fusion_file)) {
         die "unable to copy $tumor_filtered_annotated_fusion_file";
      }
  }
}

sub clonality_dir {
  my ($self, $build) = @_;

  my $patient_dir = $self->patient_dir($build);
  my $clonality_dir = $patient_dir . "/clonality";

  unless(-d $clonality_dir) {
      die $self->error_message("Clonality directory does not exist. Expected: $clonality_dir");
  }
  return $clonality_dir;
}

sub varscan_formatted_readcount_file {
  my ($self, $build) = @_;
  my $clonality_dir = $self->clonality_dir($build);
  my $readcount_file = $clonality_dir ."/allsnvs.hq.novel.tier123.v2.bed.adapted.readcounts.varscan";
  unless(-e $readcount_file) {
      die $self->error_message("Unable to find varscan formatted readcount file. Expected: $readcount_file");
  }
  return $readcount_file;
}

sub cnaseq_hmm_file {
  my ($self, $build) = @_;
  my $clonality_dir = $self->clonality_dir($build);
  my $hmm_file = $clonality_dir . "/cnaseq.cnvhmm";
  unless(-e $hmm_file) {
      die $self->error_message("Unable to find cnaseq hmm file. Expected: $hmm_file");
  }
  return $hmm_file;
}

sub has_microarray_build {
  my $self = shift;
  my $base_model;
  if ($self->exome_model) {
      $base_model = $self->exome_model;
  } elsif ($self->wgs_model) {
      $base_model = $self->wgs_model;
  } else {
      return 0;
  }
  if ($base_model->class eq 'Genome::Model::SomaticValidation') { return 0; }
  if ($base_model->tumor_model->genotype_microarray && $base_model->normal_model->genotype_microarray) {
    if($base_model->tumor_model->genotype_microarray->last_succeeded_build && $base_model->normal_model->genotype_microarray->last_succeeded_build) {
      return 1;
    } else {
      return 0;
    }
  } else {
      return 0;
  }
}

sub _get_docm_variants_file {
  my $self = shift;
  my $cancer_annotation_db = shift;
  my $cad_data_directory = $cancer_annotation_db->data_directory;
  my $docm_variants_file = $cad_data_directory . "/DOCM/DOCM_v0.1.tsv";
  if(-e $docm_variants_file ) {
    return $docm_variants_file;
  } else {
    $self->status_message("Unable to find DOCM variants file in cancer annotation db directory $docm_variants_file");
    return 0;
  }
}

1;
