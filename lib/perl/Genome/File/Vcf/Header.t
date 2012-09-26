#!/usr/bin/env perl

use above 'Genome';
use Data::Dumper;
use Test::More;

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
}

my $pkg = "Genome::File::Vcf::Header";

use_ok($pkg);

my $header_txt = <<EOS;
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="Passed all filters">
##INFO=<ID=CALLER,Number=.,Type=String,Description="Variant caller">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
EOS
my @lines = split("\n", $header_txt);
my $header = $pkg->create(lines => \@lines);

is($header->fileformat, "VCFv4.1", "fileformat");
is_deeply([$header->sample_names], ["S1", "S2", "S3"], "Sample names parsed");
is_deeply($header->info_types->{CALLER},
    {
        id => "CALLER",
        number => ".",
        type => "String",
        description => "Variant caller",
    }
    , "caller info field");

is_deeply($header->format_types->{GT},
    {
        id => "GT",
        number => "1",
        type => "String",
        description => "Genotype",
    }
    , "caller info field");

done_testing();
