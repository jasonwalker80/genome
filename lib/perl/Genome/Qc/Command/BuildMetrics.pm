package Genome::Qc::Command::BuildMetrics;

use strict;
use warnings;

use Genome;
use Data::Dumper;

class Genome::Qc::Command::BuildMetrics {
    is => 'Command::V2',
    has => [
        builds => {
            is => 'Genome::Model::Build',
            is_many => 1,
            shell_args_position => 1,
            doc => 'The builds to report QC metrics for.',
        },
        merged_metrics_file => {
            is => 'Text',
            doc => 'The file path to output build merged QC metrics.',
        },
        read_group_metrics_file => {
            is => 'Text',
            doc => 'The file path to output build read group QC metrics.',
            is_optional => 1,
        }
    ],
};
sub help_brief {
    'A command to print the QC result metrics for input builds.'
}

sub help_synopsis {
    'Dump QC result metrics from the database to tab-delimited output.'
}

sub help_detail{
    return <<"EOS"
The QC framework stores result metrics in the database for each QC result.  This tool will dump the QC result metrics for all input builds.  A tab-delimited output of all QC metrics along with build id and instrument data ids are output to the terminal.
EOS
}

sub execute {
    my $self = shift;

    my @merged_metrics;
    my @read_group_metrics;

    for my $build ($self->builds) {
        my ($merged_metrics,$read_group_metrics) = $self->metrics_for_build($build);
        push @merged_metrics, @{$merged_metrics};
        push @read_group_metrics, @{$read_group_metrics};
    }

    unless (@merged_metrics) {
        $self->error_message('Failed to find merged QC results for builds!');
        die($self->error_message);
    }

    my @merged_metrics_headers = List::MoreUtils::uniq(sort map { keys %$_ } @merged_metrics);

    my $merged_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        separator => "\t",
        headers => \@merged_metrics_headers,
        output => $self->merged_metrics_file,
    );

    for (@merged_metrics) {
        $merged_writer->write_one($_);
    }

    $merged_writer->output->close;
    if ($self->read_group_metrics_file) {
        unless (@read_group_metrics) {
            $self->error_message('Failed to find read group QC results for builds!');
            die($self->error_message);
        }
        my @read_group_metrics_headers = List::MoreUtils::uniq(sort map { keys %$_ } @read_group_metrics);

        my $read_group_writer = Genome::Utility::IO::SeparatedValueWriter->create(
            separator => "\t",
            headers => \@read_group_metrics_headers,
            output => $self->read_group_metrics_file,
        );
        for (@read_group_metrics) {
            $read_group_writer->write_one($_);
        }
        $read_group_writer->output->close;
    }

    return 1;
}

sub metrics_for_build {
    my $self = shift;
    my $build = shift;

    my @merged_metrics;
    my @read_group_metrics;
    my @qc_results = grep {$_->isa('Genome::Qc::Result')} $build->results;
    for my $qc_result (@qc_results) {
        my $as = $qc_result->alignment_result;
        my @instrument_data = $as->instrument_data;
        my %instrument_data_ids = map {$_->id => 1} @instrument_data;

        my %metrics = $qc_result->get_metrics;
        my %all_read_group_metrics;
        for my $key (keys %metrics) {
            my @metric_key_parts = split('-',$key);
            if (@metric_key_parts > 1) {
                my $metric_label = $metric_key_parts[0];
                my $metric_name = $metric_key_parts[1];
                if ($instrument_data_ids{$metric_label}) {
                    my $metric_value = delete($metrics{$key});
                    $all_read_group_metrics{$metric_label}{$metric_name} = $metric_value;
                }
            }
        }
        for my $read_group_id (keys %all_read_group_metrics) {
            my %read_group_metrics = %{$all_read_group_metrics{$read_group_id}};
            $read_group_metrics{read_group_id} = $read_group_id;
            $read_group_metrics{build_id} = $build->id;
            push @read_group_metrics, \%read_group_metrics;
        }
        $metrics{build_id} = $build->id;
        $metrics{read_groups} = scalar(keys %instrument_data_ids);
        $metrics{read_group_ids} = join(',',keys %instrument_data_ids);
        push @merged_metrics, \%metrics;
    }

    return (\@merged_metrics, \@read_group_metrics);
}


1;
