#!/usr/bin/env genome-perl

BEGIN { 
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

use strict;
use warnings;

use above "Genome";
use Test::Exception;
use Test::More;
use Genome::VariantReporting::BamReadcount::TestHelper qw(bam_readcount_line
    create_entry create_deletion_entry bam_readcount_line_deletion);

my $pkg = "Genome::VariantReporting::BamReadcount::VafFilter";
use_ok($pkg);
my $factory = Genome::VariantReporting::Framework::Factory->create();
isa_ok($factory->get_class('filters', $pkg->name), $pkg);

subtest "__errors__ missing parameters" => sub {
    my $filter = $pkg->create(sample_name => "S1");
    my @errors = $filter->__errors__;
    is(scalar(@errors), 1, "One error found");
    like($errors[0]->desc, qr/Must define at least one of min_vaf or max_vaf/, "Error message as expected");
};

subtest "__errors__ invalid parameters" => sub {
    my $filter = $pkg->create(min_vaf => 100, max_vaf => 90, sample_name => "S1");
    my @errors = $filter->__errors__;
    is(scalar(@errors), 1, "One error found");
    like($errors[0]->desc, qr/Max_vaf must be larger or equal to min_vaf/, "Error message as expected");
};

subtest "pass min vaf" => sub {
    my $min_vaf = 90;
    my $filter = $pkg->create(min_vaf => $min_vaf, sample_name => "S1");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 1,
        C => 0,
        AA => 0,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry passes filter with min_vaf $min_vaf");
};

subtest "min vaf for insertion" => sub {
    my $min_vaf = 5;
    my $filter = $pkg->create(min_vaf => $min_vaf, sample_name => "S4");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 0,
        C => 0,
        AA => 1,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry passes filter with min_vaf $min_vaf");
};

subtest "min vaf for deletion" => sub {
    my $min_vaf = 5;
    my $filter = $pkg->create(min_vaf => $min_vaf, sample_name => "S1");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        A => 1,
    );
    my $entry = create_deletion_entry(bam_readcount_line_deletion);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry passes filter with min_vaf $min_vaf");
};

subtest "pass max vaf" => sub {
    my $max_vaf = 100;
    my $filter = $pkg->create(max_vaf => $max_vaf, sample_name => "S2");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 1,
        C => 1,
        AA => 0,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry passes filter with max_vaf $max_vaf");
};

subtest "Bamreadcount entry is a ." => sub {
    my $filter = $pkg->create(max_vaf => 100, sample_name => "S2");
    my %expected_return_values = (
        G => 1,
        C => 1,
        AA => 0,
    );
    my $entry = create_entry(".");
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry passes filter when bam-readcount not available");
};

subtest "pass min and max vaf" => sub {
    my $min_vaf = 90;
    my $max_vaf = 100;
    my $filter = $pkg->create(min_vaf => $min_vaf, max_vaf => $max_vaf, sample_name => "S2");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 1,
        C => 0,
        AA => 0,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry passes filter with min_vaf $min_vaf and max_vaf $max_vaf");
};

subtest "fail min vaf" => sub {
    my $min_vaf = 100;
    my $filter = $pkg->create(min_vaf => $min_vaf, sample_name => "S1");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 0,
        C => 0,
        AA => 0,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry fails filter with min_vaf $min_vaf");
};

subtest "fail max vaf" => sub {
    my $max_vaf = 90;
    my $filter = $pkg->create(max_vaf => $max_vaf, sample_name => "S2");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 0,
        C => 1,
        AA => 0,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry fails filter with max_vaf $max_vaf");
};

subtest "fail heterozygous non-reference sample" => sub {
    my $min_vaf = 90;
    my $filter = $pkg->create(min_vaf => $min_vaf, sample_name => "S2");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 1,
        C => 0,
        AA => 0,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry fails filter with min_vaf $min_vaf");
    cmp_ok(Genome::VariantReporting::BamReadcount::VafCalculator::calculate_vaf(
        $filter->get_readcount_entry($entry), 'C', 'A'), '<', 0.3, "VAF is very low");
    cmp_ok(Genome::VariantReporting::BamReadcount::VafCalculator::calculate_vaf(
        $filter->get_readcount_entry($entry), 'G', 'A'), '>', 90, "VAF is high");
};

subtest "pass heterozygous non-reference sample" => sub {
    my $min_vaf = 0.02;
    my $filter = $pkg->create(min_vaf => $min_vaf, sample_name => "S2");
    lives_ok(sub {$filter->validate}, "Filter validates");

    my %expected_return_values = (
        G => 1,
        C => 1,
        AA => 0,
    );
    my $entry = create_entry(bam_readcount_line);
    is_deeply({$filter->filter_entry($entry)}, \%expected_return_values, "Entry passes filter with min_vaf $min_vaf");
};

subtest "no bam readcount entry" => sub {
    my $no_readcount_entry = create_entry();
    my $filter = $pkg->create(
        sample_name => 'S1',
        min_vaf => 5,
    );
    lives_ok(sub {$filter->validate}, "Filter validates ok");
    my %expected_return_values = (
        G => 1,
        C => 0,
        AA => 0,
    );
    is_deeply({$filter->filter_entry($no_readcount_entry)}, \%expected_return_values, "Sample 1 return values as expected");
};

done_testing;
