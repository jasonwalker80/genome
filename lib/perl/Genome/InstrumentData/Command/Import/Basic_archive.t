#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above "Genome";

require Genome::Utility::Test;
require File::Compare;
use Test::More;

use_ok('Genome::InstrumentData::Command::Import::Basic') or die;
use_ok('Genome::InstrumentData::Command::Import::WorkFlow::Helpers') or die;

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'fastq/v1');
my $source_file = $test_dir.'/input.fastq.tgz';
ok($source_file, 'source archive exists');

my $analysis_project = Genome::Config::AnalysisProject->create(name => '__TEST_AP__');
ok($analysis_project, 'create analysis project');
my $library = Genome::Library->create(
    name => '__TEST_SAMPLE__-extlibs', sample => Genome::Sample->create(name => '__TEST_SAMPLE__')
);
ok($library, 'Create library');

my $cmd = Genome::InstrumentData::Command::Import::Basic->create(
    analysis_project => $analysis_project,
    library => $library,
    source_files => [$source_file],
    import_source_name => 'broad',
    instrument_data_properties => [qw/ lane=2 flow_cell_id=XXXXXX /],
    original_format => 'fastq',
);
ok($cmd, "create import command");
ok($cmd->execute, "excute import command");

my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($source_file.'.md5');
ok($md5, 'load source md5');
my @instrument_data = map { $_->instrument_data } Genome::InstrumentDataAttribute->get(
    attribute_label => 'original_data_path_md5',
    attribute_value => $md5,
);
is(@instrument_data, 1, "got instrument data for md5 $md5") or die;;
my $instrument_data = $instrument_data[0];
is($instrument_data->original_data_path, $source_file, 'original_data_path correctly set');
is($instrument_data->import_format, 'bam', 'import_format is bam');
is($instrument_data->sequencing_platform, 'solexa', 'sequencing_platform correctly set');
is($instrument_data->is_paired_end, 1, 'is_paired_end correctly set');
is($instrument_data->read_count, 2000, 'read_count correctly set');
is(eval{ $instrument_data->attributes(attribute_label => 'lane')->attribute_value }, 2, 'lane correctly set');
is(eval{ $instrument_data->attributes(attribute_label => 'flow_cell_id')->attribute_value }, 'XXXXXX', 'flow_cell_id correctly set');
is($instrument_data->analysis_projects, $analysis_project, 'set analysis project');

my $bam_path = $instrument_data->bam_path;
ok(-s $bam_path, 'bam path exists');
is($bam_path, $instrument_data->data_directory.'/all_sequences.bam', 'bam path correctly named');
is(eval{$instrument_data->attributes(attribute_label => 'bam_path')->attribute_value}, $bam_path, 'set attributes bam path');
is(File::Compare::compare($bam_path, $test_dir.'/input.fastq.bam'), 0, 'bam matches');
is(File::Compare::compare($bam_path.'.flagstat', $test_dir.'/input.fastq.bam.flagstat'), 0, 'flagstat matches');

my $allocation = $instrument_data->disk_allocation;
ok($allocation, 'got allocation');
ok($allocation->kilobytes_requested > 0, 'allocation kb was set');

#print $instrument_data->data_directory."\n";<STDIN>;
done_testing();
