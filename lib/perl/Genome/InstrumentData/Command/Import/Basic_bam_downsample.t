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

my $analysis_project = Genome::Config::AnalysisProject->create(name => '__TEST_AP__');
ok($analysis_project, 'create analysis project');

my $library = Genome::Library->create(
    name => '__TEST_SAMPLE__-extlibs', sample => Genome::Sample->create(name => '__TEST_SAMPLE__')
);
ok($library, 'Create library');

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam/v4');
my $source_bam = $test_dir.'/test.bam';
ok(-s $source_bam, 'source bam exists') or die;

my $cmd = Genome::InstrumentData::Command::Import::Basic->create(
    analysis_project => $analysis_project,
    library => $library,
    source_files => [$source_bam],
    import_source_name => 'broad',
    downsample_ratio => .25,
    instrument_data_properties => [qw/ lane=2 flow_cell_id=XXXXXX /],
);
ok($cmd, "create import command");
ok($cmd->execute, "excute import command");

my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($source_bam.'.md5');
ok($md5, 'load source md5');
my @instrument_data = map { $_->instrument_data } Genome::InstrumentDataAttribute->get(
    attribute_label => 'original_data_path_md5',
    attribute_value => $md5,
);
is(@instrument_data, 1, "got instrument data for md5 $md5") or die;;
my $instrument_data = $instrument_data[0];
is($instrument_data->original_data_path, $source_bam, 'original_data_path correctly set');
is($instrument_data->import_format, 'bam', 'import_format is bam');
is($instrument_data->sequencing_platform, 'solexa', 'sequencing_platform correctly set');
is($instrument_data->is_paired_end, 1, 'is_paired_end correctly set');
is($instrument_data->read_count, 68, 'read_count correctly set');
is($instrument_data->read_length, 100, 'read_length correctly set');
cmp_ok(eval{$instrument_data->attributes(attribute_label => 'downsample_ratio')->attribute_value;}, '==', 0.25, 'downsample_ratio correctly set');
is(eval{$instrument_data->attributes(attribute_label => 'segment_id')->attribute_value;}, 2883581797, 'segment_id correctly set');
is(eval{$instrument_data->attributes(attribute_label => 'original_data_path_md5')->attribute_value;}, 'f81fbc3d3a6b57d11e60b016bb2c950c', 'original_data_path_md5 correctly set');
is($instrument_data->analysis_projects, $analysis_project, 'set analysis project');

my $bam_path = $instrument_data->bam_path;
ok(-s $bam_path, 'bam path exists');
is($bam_path, File::Spec->catfile($instrument_data->data_directory, 'all_sequences.bam'), 'bam path correctly named');
is(eval{$instrument_data->attributes(attribute_label => 'bam_path')->attribute_value}, $bam_path, 'set attributes bam path');
is(# compare to the split version b/c just the downsampled bam header has attrs sorted differently
    File::Compare::compare($bam_path, File::Spec->catfile($test_dir, 'test.clean.sorted.downsampled.split-by-rg.bam')),
    0, 'bam matches',
);
is(# compare to the downsampled flagstat
    File::Compare::compare($bam_path.'.flagstat', File::Spec->catfile($test_dir, 'test.clean.sorted.downsampled.bam.flagstat')),
    0, 'flagstat matches',
);

my $allocation = $instrument_data->disk_allocation;
ok($allocation, 'got allocation');
ok($allocation->kilobytes_requested > 0, 'allocation kb was set');

done_testing();