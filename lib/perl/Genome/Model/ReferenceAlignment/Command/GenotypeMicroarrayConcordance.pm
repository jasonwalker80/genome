package Genome::Model::ReferenceAlignment::Command::GenotypeMicroarrayConcordance;

use strict;
use warnings;

use Genome;

use Params::Validate ':types';
use File::Spec;

class Genome::Model::ReferenceAlignment::Command::GenotypeMicroarrayConcordance {
    is => 'Command::V2',
    doc => 'Determine concordance between a genotype microarray VCF and a Reference Alignment build',
    has => [
        build => {
            is => 'Genome::Model::Build::ReferenceAlignment',
            doc => 'use this build to find genotype microarray builds and report concordance for all instrument data',
            shell_args_position => 1,
        },
        minimum_depth => {
            is => 'Integer',
            doc => 'The minimum depth required to consider concordant.',
            example_values => ['4'],
        },
        bedtools_version => {
            is => 'Text',
            doc => 'Version of bed tools intersect to use. Recommended version: '.Genome::Model::Tools::BedTools->latest_bedtools_version,
            valid_values => [ Genome::Model::Tools::BedTools->available_bedtools_versions ],
        },
        picard_version => {
            is => 'Text',
            doc => 'The version of picard to use.',
            valid_values => [ Genome::Model::Tools::Picard::GenotypeConcordance->available_picard_versions ],
        },
        dbsnp_build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'The dbsnp build used to generate the microarray VCFs.',
            example_values => ['141765601'],
        },
    ],
    has_optional =>[
        genotype_microarray_sample => {
            is => 'Genome::Sample',
            doc => 'The sample for which genotype microarray data exists to compare to the sequenced sample.',
        },
        output_dir => {
            is => 'Text',
            doc => 'A directory to output the raw Picard files to.',
        },
        _seqdict => {},
    ],
};

sub help_detail {
    return "Generate picard genotype concordance metrics between the genotype microarray VCF and a reference alignment build.";
}

sub execute {
    my $self = shift;

    # The sequence dictionary is required by Picard
    $self->resolve_sequence_dictionary;

    # The sequence data SNVs VCF from the RefAlign Build
    my $build_vcf = $self->resolve_snvs_vcf_for_build();

    # Determine the sample to compare genotype microarray data with
    $self->resolve_genotype_microarray_sample;

    # The microarray genotypes as a VCF file
    my $microarray_vcf = $self->resolve_genotype_microarray_vcf;

    # Run Picard GenotypeConcordance and return the output prefix
    my $output_prefix = $self->run_picard_genotype_concordance($build_vcf,$microarray_vcf);

    # Parse the output metrics into a hash ref
    my $data = $self->output_summary_hash_ref($output_prefix);

    # Write the summary Picard GC metrics along with some build summary metrics
    $self->write_output_summary($data);

    return 1;
}

sub picard_metrics {
    my $class = shift;
    return qw/
                 HET_SENSITIVITY
                 HET_PPV
                 HOMVAR_SENSITIVITY
                 HOMVAR_PPV
                 VAR_SENSITIVITY
                 VAR_PPV
             /;
}

sub summary_headers {
    my $class = shift;
    my @headers = qw/
                        instrument_data_id
                        flow_cell_id
                        lane
                        clusters
                        sample_name
                        microarray_sample_name
                        snvs_build_id
                    /;
    push @headers, $class->picard_metrics;
    return @headers;
}

sub resolve_snvs_vcf_for_build {
    my $self = shift;

    my $build = $self->build;
    my $build_vcf = $build->get_detailed_snvs_vcf;
    unless (-e $build_vcf) {
        $self->error_message('Unable to find the snvs VCF for build : '. $build->__display_name__);
        die($self->error_message);
    }
    return $build_vcf;
}

sub resolve_genotype_microarray_vcf {
    my $self = shift;

    # Get the sorted VCF file name
    my $sorted_microarray_vcf = $self->sorted_microarray_vcf_for_genotype_microarray_sample();

    $self->debug_message('Get or create genotype VCF result for sample: '. $self->genotype_microarray_sample->__display_name__);

    my $vcf_result = Genome::InstrumentData::Microarray::Result::Vcf->get_or_create(
        sample => $self->genotype_microarray_sample,
        known_sites_build => $self->dbsnp_build,
        users => Genome::SoftwareResult::User->user_hash_for_build($self->build),
    );
    my $microarray_vcf = $vcf_result->vcf_path;
    unless (-e $microarray_vcf) {
        $self->error_message('Failed to get or create microarray VCF for sample: '. $self->genotype_microarray_sample->__display_name__);
        die($self->error_message);
    }
    $self->sort_vcf($microarray_vcf,$sorted_microarray_vcf);
    return $sorted_microarray_vcf;
}

sub sorted_microarray_vcf_for_genotype_microarray_sample {
    my $self = shift;

    my $genotype_microarray_sample = $self->genotype_microarray_sample;
    # Create the sorted VCF file name for this sample
    my $sorted_microarray_vcf = $self->create_temp_file_path(
        Genome::Utility::Text::sanitize_string_for_filesystem($genotype_microarray_sample->name).'_microarray_sorted.vcf.gz'
      );

    return $sorted_microarray_vcf;
}

sub sort_vcf {
    my $self = shift;
    my $input_vcf = shift;
    my $output_vcf = shift;

    # Sort the VCF which adds the contig info to the VCF header (required by Picard GenotypeConcordance)
    # A suffix of .vcf.gz will also compress and index the resulting VCF

    my $sort_cmd = Genome::Model::Tools::Picard::SortVcf->create(
        input_vcf => $input_vcf,
        output_vcf => $output_vcf,
        sequence_dictionary => $self->_seqdict,
        use_version => $self->picard_version,
    );
    unless ($sort_cmd) {
        $self->error_message('Failed to create sort VCF command!');
        die($self->error_message);
    }
    unless ($sort_cmd->execute) {
        $self->error_message('Failed to sort VCF: '. $input_vcf);
        die($self->error_message);
    }
    return 1;
}

sub create_temp_file_path {
    my $self = shift;
    my $basename = shift;

    if ($self->output_dir) {
        return File::Spec->join($self->output_dir,$basename);
    } else {
        return Genome::Sys->create_temp_file_path($basename);
    }
}

sub intersect_vcfs {
    my $self = shift;
    my ($a_vcf,$b_vcf,$basename) = @_;
    
    my $raw_vcf = $self->create_temp_file_path($basename .'.vcf');
    my $intersect_cmd = Genome::Model::Tools::BedTools::Intersect->create(
        input_file_a => $a_vcf,
        input_file_a_format => 'vcf',
        input_file_b => $b_vcf,
        output_file => $raw_vcf,
        header => 1,
        use_version => $self->bedtools_version,
    );
    unless ($intersect_cmd) {
        $self->error_message('Failed to create bedtools intersect!');
        die($self->error_message);
    }
    unless ($intersect_cmd->execute) {
        $self->error_message('Failed to execute bedtools intersect!');
        die($self->error_message);
    }
    
    my $sorted_vcf = $self->create_temp_file_path($basename .'.vcf.gz');
    $self->sort_vcf($raw_vcf,$sorted_vcf);

    return $sorted_vcf;
}

sub resolve_roi_intervals_from_build {
    my $self = shift;
    my $build = $self->build;
    
    my $roi_bed = $self->create_temp_file_path($build->id .'.bed');
    $build->region_of_interest_set_bed_file($roi_bed);

    my $roi_intervals = $self->create_temp_file_path($build->id .'.intervals');
    my $bed_to_intervals_cmd = Genome::Model::Tools::Picard::BedToIntervalList->create(
        input => $roi_bed,
        output => $roi_intervals,
        sequence_dictionary => $self->_seqdict,
        use_version => $self->picard_version,
    );
    unless ($bed_to_intervals_cmd) {
        $self->error_message('Failed to create BedToIntervals!');
        die($self->error_message);
    }
    unless ($bed_to_intervals_cmd->execute) {
        $self->error_message('Failed to execute BedToIntervals!');
        die($self->error_message);
    }
    return $roi_intervals;
}

sub resolve_genotype_microarray_sample {
    my $self = shift;

    return 1 if defined($self->genotype_microarray_sample);
    $self->debug_message('Attempting to resolve the genotype sample from the input RefAlign build: '. $self->build->id);

    my $microarray_build = $self->build->genotype_microarray_build;
    if ($microarray_build) {
        my $genotype_microarray_sample = $microarray_build->model->subject;
        if (defined($genotype_microarray_sample)) {
            $self->debug_message('Resolved genotype microarray sample from microarray build: '. $genotype_microarray_sample->__display_name__);
            $self->genotype_microarray_sample($genotype_microarray_sample);
        } else {
            $self->error_message('Failed to find the model subject for microarray build '. $microarray_build->id);
        }
    } else {
        my $genotype_microarray_sample = $self->build->model->subject;
        if (defined($genotype_microarray_sample)) {
            $self->debug_message('Resolved genotype microarray sample from RefAlign build: '. $genotype_microarray_sample->__display_name__);
            $self->genotype_microarray_sample($genotype_microarray_sample);
        } else {
            $self->error_message('Failed to find the model subject for RefAlign build '. $self->build->id);
        }
    }
    return 1;
}

sub resolve_sequence_dictionary {
    my $self = shift;

    my $ref_seq_build = $self->build->reference_sequence_build;
    my $seqdict = $ref_seq_build->get_sequence_dictionary('sam',$ref_seq_build->species_name,$self->picard_version);
    $self->_seqdict($seqdict);

    return 1;
}

sub run_picard_genotype_concordance {
    my $self = shift;
    my $build_vcf = shift;
    my $microarray_vcf = shift;
    
    # Without an intersection the PPV is really high.  Presumably caused by many False-Positive SNVs
    my $intersect_vcf_basename = Genome::Utility::Text::sanitize_string_for_filesystem($self->genotype_microarray_sample->name) .'_x_'. $self->build->id;
    my $intersect_vcf = $self->intersect_vcfs($build_vcf,$microarray_vcf,$intersect_vcf_basename);

    my $output_prefix = $self->create_temp_file_path($self->build->id);
    my %gc_params = (
        truth_vcf => $microarray_vcf,
        call_vcf => $intersect_vcf,
        output => $output_prefix,
        truth_sample => $self->genotype_microarray_sample->name,
        call_sample => $self->build->model->subject->name_in_vcf,
        min_dp => $self->minimum_depth,
        use_version => $self->picard_version,
    );
    # For exome limit to ROI from input build
    if ($self->build->model->is_capture) {
        $gc_params{intervals} = $self->resolve_roi_intervals_from_build();
    }
    my $gc_cmd = Genome::Model::Tools::Picard::GenotypeConcordance->create(%gc_params);
    unless ($gc_cmd) {
        $self->error_message('Failed to create GenotypeConcordance!');
        die($self->error_message);
    }
    unless ($gc_cmd->execute) {
        $self->error_message('Failed to execute GenotypeConcordance!');
        die($self->error_message);
    }
    return $output_prefix;
}

sub output_summary_hash_ref {
    my $self = shift;
    my $output_prefix = shift;

    
    my @instrument_data = $self->build->instrument_data;
    my $instrument_data_ids = join(',', map {$_->id} @instrument_data);
    my $flow_cell_ids = join(',', map {$_->flow_cell_id} @instrument_data);
    my $lanes = join(',', map {$_->lane} @instrument_data);
    my $clusters = join(',', map {$_->clusters} @instrument_data);
    my %data = (
        snvs_build_id => $self->build->id,
        instrument_data_id => $instrument_data_ids,
        sample_name => $self->build->model->subject->name,
        microarray_sample_name => $self->genotype_microarray_sample->name,
        flow_cell_id => $flow_cell_ids,
        lane => $lanes,
        clusters => $clusters,
    );
    
    my $gc_summary_metrics_file = $output_prefix .'.genotype_concordance_summary_metrics';
    if (-e $gc_summary_metrics_file) {
        my $gc_summary_metrics_hash_ref = Genome::Model::Tools::Picard->parse_file_into_metrics_hashref($gc_summary_metrics_file, 'VARIANT_TYPE');
        my $gc_summary_snp_metrics = $gc_summary_metrics_hash_ref->{'SNP'};
        for my $key ($self->picard_metrics) {
            $data{$key} = $gc_summary_snp_metrics->{$key};
        }
    } else {
        die('Failed to find GenotypeConcordance summary file: '. $gc_summary_metrics_file);
    }

    return \%data;
}

sub write_output_summary {
    my $self = shift;
    my $data = shift;
    
    my @headers = $self->summary_headers;
    my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@headers,
        in_place_of_null_value => 'na',
    );
    $writer->write_one($data);
    $writer->output->close;
}

1;

