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

my $test_dir = Genome::Utility::Test->data_dir_ok('Genome::InstrumentData::Command::Import', 'bam-rg-multi/v3');
my $source_bam = $test_dir.'/input.rg-multi.bam';
ok(-s $source_bam, 'source bam exists') or die;

my $cmd = Genome::InstrumentData::Command::Import::Basic->create(
    analysis_project => $analysis_project,
    library => $library,
    source_files => [$source_bam],
    import_source_name => 'broad',
    instrument_data_properties => [qw/ lane=2 flow_cell_id=XXXXXX /],
);
ok($cmd, "create import command");
ok($cmd->execute, "excute import command");

my $md5 = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->load_md5($source_bam.'.md5');
ok($md5, 'load source md5');
my @instrument_data_attributes = Genome::InstrumentDataAttribute->get(
    attribute_label => 'original_data_path_md5',
    attribute_value => $md5,
);
my @instrument_data = Genome::InstrumentData::Imported->get(id => [ map { $_->instrument_data_id } @instrument_data_attributes ]);
is(@instrument_data, 4, "got instrument data for md5 $md5");

my %read_groups = (
    2883581797 => { paired => 34, singleton => 94 },
    2883581798 => { paired => 36, singleton => 92 },
);

for my $instrument_data ( @instrument_data ) {
    my $read_group_id = $instrument_data->attributes(attribute_label => 'segment_id')->attribute_value;

    my $paired_key = $instrument_data->is_paired_end? 'paired' : 'singleton';
    ok(!$read_groups{$read_group_id}{$paired_key.'-seen'}, "read group $read_group_id $paired_key not seen");
    $read_groups{$read_group_id}{$paired_key.'-seen'} = 1;

    is($instrument_data->original_data_path, $source_bam, 'original_data_path correctly set');
    is($instrument_data->import_format, 'bam', 'import_format is bam');
    is($instrument_data->sequencing_platform, 'solexa', 'sequencing_platform correctly set');
    is($instrument_data->read_count, $read_groups{$read_group_id}{$paired_key}, 'read_count correctly set');
    is($instrument_data->read_length, 100, 'read_length correctly set');
    is($instrument_data->analysis_projects, $analysis_project, 'set analysis project');

    my $bam_path = $instrument_data->bam_path;
    ok(-s $bam_path, 'bam path exists');
    is($bam_path, $instrument_data->data_directory.'/all_sequences.bam', 'bam path correctly named');
    is(eval{$instrument_data->attributes(attribute_label => 'bam_path')->attribute_value}, $bam_path, 'set attributes bam path');
    is(File::Compare::compare($bam_path, $test_dir.'/'.join('.', $read_group_id, $paired_key,'bam')), 0, 'bam matches');
    is(File::Compare::compare($bam_path.'.flagstat', $test_dir.'/'.join('.', $read_group_id, $paired_key, 'bam.flagstat')), 0, 'flagstat matches');

    my $allocation = $instrument_data->disk_allocation;
    ok($allocation, 'got allocation');
    ok($allocation->kilobytes_requested > 0, 'allocation kb was set');

    #print join("\t", $read_group_id, $paired_key, $instrument_data->bam_path), "\n";
}

done_testing();
