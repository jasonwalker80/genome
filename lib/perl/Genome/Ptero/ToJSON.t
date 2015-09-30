#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => 'all';

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

use Test::More;
use above 'Genome';
use Genome::Utility::Test qw(compare_ok);
use File::Basename qw(dirname basename);
use File::Spec qw();
use File::Copy qw(copy);


for my $test_directory (glob test_data_directory('*')) {
    my $test_name = basename($test_directory);
    my $workflow = Genome::WorkflowBuilder::DAG->from_xml_filename(workflow_xml_file($test_name));

    my $json_filename = do {
        local %ENV; # otherwise the JSON contains the current user's environment.
        $ENV{XGENOME_PTERO_SHELL_COMMAND_SERVICE_URL} = 'http://example.com/v1';
        $ENV{XGENOME_PTERO_LSF_SERVICE_URL} = 'http://lsf.example.com/v1';
        my $ptero_builder = $workflow->get_ptero_builder();

        my $json_filename = Genome::Sys->create_temp_file_path();
        Genome::Sys->write_file($json_filename, $ptero_builder->to_json());

        $json_filename;
    };

    my $expected_json_filename = expected_workflow_json_file($test_name);
    if ($ENV{GENERATE_TEST_DATA}) {
        copy($json_filename, $expected_json_filename);
    }

    compare_ok($json_filename, $expected_json_filename,
        "$test_name JSON looks ok",
        filters => [
            qr/^.*"workingDirectory" :.*$/,
            qr/^.*"user" :.*$/,
            qr/^.*"cwd" :.*$/,
        ],
    );
}


done_testing();


sub workflow_xml_file {
    my $name = shift;

    return File::Spec->join(test_data_directory($name), 'workflow.xml');
}

sub expected_workflow_json_file {
    my $name = shift;

    return File::Spec->join(test_data_directory($name), 'ptero_workflow.json');
}

sub test_data_directory {
    my $name = shift;

    return File::Spec->join(dirname(__FILE__), 'workflow_tests', $name);
}
