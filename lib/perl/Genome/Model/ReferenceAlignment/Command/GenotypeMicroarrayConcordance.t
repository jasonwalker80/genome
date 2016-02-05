#!/usr/bin/env genome-perl

use strict;
use warnings 'FATAL';

use above 'Genome';

use Genome::Test::Factory::Sample;
use Test::More;
plan tests => 3;

my $class = 'Genome::Model::ReferenceAlignment::Command::GenotypeMicroarrayConcordance';
use_ok($class) or die;

my $sample = Genome::Test::Factory::Sample->generate_obj(name => '_TEST(1)_SAMPLE_');

my $gc = $class->create(
   genotype_microarray_sample => $sample,
);
ok($gc, 'create genotype concordance command');

my $sorted_vcf = $gc->sorted_microarray_vcf_for_genotype_microarray_sample();
my $sanitized_sample_name = Genome::Utility::Text::sanitize_string_for_filesystem($sample->name);
like($sorted_vcf, qr/$sanitized_sample_name/, 'sorted_microarray_vcf_for_genotype_microarray_sample has sanitized sample name');

done_testing();
