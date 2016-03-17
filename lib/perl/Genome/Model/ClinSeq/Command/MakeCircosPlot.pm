package Genome::Model::ClinSeq::Command::MakeCircosPlot;
use strict;
use warnings;
use Switch;
use Genome; 

# Written by Ben Ainscough and Scott Smith, based on prototype from Obi Griffith
# See JIRA issue https://jira.gsc.wustl.edu/browse/TD-691
# Modified to support exome only data by Kilannin Krysiak

class Genome::Model::ClinSeq::Command::MakeCircosPlot {
   is => 'Command::V2',
    has_input => [
        
        output_directory    => { is => 'FilesystemPath',
                                doc => 'Directory where output will be written', },
                                
        build               => { is => 'Genome::Model::Build::ClinSeq',
                                doc => 'Clinseq build' },
                                
        candidate_fusion_infile      => {  is => 'FilesystemPath', doc => 'fusion_output_file file from SummarizeSvs.pm' , is_optional => 1 },
        cnv_hmm_file				 => {  is => 'FilesystemPath', doc => 'cnv_hmm_file is from GenerateClonalityPlots.pm', is_optional => 1 },
        coding_hq_de_file			 => {  is => 'FilesystemPath', doc => 'coding_hq_de_file is from CufflinksDifferentialExpression.pm', is_optional => 1 },
        tumor_fpkm_topnpercent_file  => {  is => 'FilesystemPath', doc => 'tumor_fpkm_topnpercent_file is from CufflinksExpressionAbsolute.pm', is_optional => 1 },
        import_snvs_indels_result    => {  is => 'Boolean', doc => 'Used in the to link in workflow', is_optional => 1 },
        annotated_variants_tsv    => {
            is => 'Text',
            doc => 'A TSV file of annotated variants.',
            is_optional => 1
        },
        gene_ampdel_file             => {  is => 'FilesystemPath', doc => 'gene_ampdel_file is from RunCnView.pm and is used to provide gene lables for Deletions and Focal Amps', is_optional => 1 },

        #TODO: Define all input files as optional inputs here
        #TODO: Each of these will have to be defined as output on the commands they come from.


    ], 
    has_param => [
        use_version         => { is => 'Text',
                                valid_values => [ Genome::Sys->sw_versions("circos") ],
                                default_value => '0.64',
                                doc => 'the version of circos to use' },
    ],
    doc => 'interfaces with the circos program to produce a circos plot for a clin-seq build',
};

sub sub_command_category { 'pipeline' }

sub help_detail {
  return <<EOS
Generate a Circos plot for a clinseq build.
EOS
}

sub help_synopsis {
  return <<EOS
    genome model clin-seq make-circos-plot --build=aee9a31051754702a9c2835d63abc812 --output-directory /tmp/outdir
EOS
}

sub rnaExpression {
	my $self = shift;
	my $build = shift;
	my $output_directory = shift;
	my $dataDir = shift;
	my $genes_AmpDel = shift;
	my $genes_noAmpDel = shift;
# 	my $config = shift;
	my $RNAconfig = "";
	
	###Differential Expression
	# Differential Expression data is only included if rnaseq builds are present for tumor and normal
	# If not the rna expression is displayed.
	if($build->normal_rnaseq_build || $build->tumor_rnaseq_build){
		if($build->normal_rnaseq_build){
			Genome::Sys->copy_file("$dataDir/rnaseq/cufflinks_differential_expression/genes/case_vs_control.coding.hq.de.tsv", "$output_directory/raw/case_vs_control.coding.hq.de.tsv");
			if($self->coding_hq_de_file){
				system("rm -f $output_directory/raw/case_vs_control.coding.hq.de.tsv");
				Genome::Sys->copy_file($self->coding_hq_de_file , "$output_directory/raw/case_vs_control.coding.hq.de.tsv");
			}
		}else{
			Genome::Sys->copy_file("$dataDir/rnaseq/tumor/cufflinks_expression_absolute/genes/genes.fpkm.expsort.top1percent.tsv", "$output_directory/raw/genes.fpkm.expsort.top1percent.tsv");
			if($self->tumor_fpkm_topnpercent_file){
				system("rm -f $output_directory/raw/genes.fpkm.expsort.top1percent.tsv");
				Genome::Sys->copy_file($self->tumor_fpkm_topnpercent_file , "$output_directory/raw/genes.fpkm.expsort.top1percent.tsv");
			}
		}
	}
	
	if($build->normal_rnaseq_build || $build->tumor_rnaseq_build){
    
    ###Differential Expression
    # Differential Expression data is only included if rnaseq builds are present for tumor and normal
    # If not the rna expression is displayed in this track.

    if($build->normal_rnaseq_build){
        my $diffExpression = Genome::Sys->read_file("$output_directory/raw/case_vs_control.coding.hq.de.tsv");
        my $diffExpressionPositive_fh = Genome::Sys->open_file_for_writing("$output_directory/data/case_vs_control.coding.hq.de.positive.txt");
        my $diffExpressionNegative_fh = Genome::Sys->open_file_for_writing("$output_directory/data/case_vs_control.coding.hq.de.negative.txt");
        while ($diffExpression =~ /ENS\w+\s(\w+)\s\S+\s\S+\s(\w+):(\d+)-(\d+)\s\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s\S+\s(\S+)/g) {
            $$genes_noAmpDel{$1}="hs$2\t$3\t$4";
			$$genes_AmpDel{$1}="hs$2\t$3\t$4";
            if($5>=2){
                print $diffExpressionPositive_fh ("hs$2 $3 $4 $5\n");
            }elsif($5<=-2){
                print $diffExpressionNegative_fh ("hs$2 $3 $4 $5\n");
            }
        }
        $diffExpressionPositive_fh->close;    
        $diffExpressionNegative_fh->close;    
        
        $RNAconfig.=<<EOS;

#DIFFERENTIAL EXPRESSION DATA 
<plot>
# The type sets the format of the track.
type = histogram
file = $output_directory/data/case_vs_control.coding.hq.de.positive.txt
min=-10
max=10 #Cap, otherwise outliers hide everything

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the ideogram.
r1   = 0.80r
r0   = 0.70r

# Histograms can have both a fill and outline. The default outline is 1px thick black.
fill_color = red

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.
thickness = 0p

# Do not join histogram bins that do not abut.
extend_bin = no
orientation = out

#Draw background to highlight track
<backgrounds>
show  = yes
<background>
color = vvlgrey
</background>
</backgrounds>

#Draw axes lines
<axes>
show = yes
thickness = 1
color     = lgrey
<axis>
spacing   = 0.1r
</axis>
</axes>

</plot>
    
<plot>
# The type sets the format of the track.
type = histogram
file = $output_directory/data/case_vs_control.coding.hq.de.negative.txt
min=-10
max=10 #Cap, otherwise outliers hide everything

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the ideogram.
r1   = 0.80r
r0   = 0.70r

# Histograms can have both a fill and outline. The default outline is 1px thick black.
fill_color = blue

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.
thickness = 0p

# Do not join histogram bins that do not abut.
extend_bin = no
orientation = out

</plot>

EOS
    
    }else{
        my $tumor_rnaseq_build=$build->tumor_rnaseq_build;
        my $reader = Genome::Utility::IO::SeparatedValueReader->create(separator => "\t", input => "$output_directory/raw/genes.fpkm.expsort.top1percent.tsv");
        my $expression_fh = Genome::Sys->open_file_for_writing("$output_directory/data/genes.fpkm.expsort.top1percent.tsv");
        while(my $data = $reader->next) {
            my ($chr, $start, $stop) = split(/[:-]/, $data->{locus});
            $$genes_noAmpDel{$data->{mapped_gene_name}}=join("\t","hs$chr", $start, $stop);
            $$genes_AmpDel{$data->{mapped_gene_name}}=join("\t","hs$chr", $start, $stop);
            #if($5>=2 || $5<=-2){ Is this necessary
                print $expression_fh ("hs$chr $start $stop ".log($data->{FPKM})/log(2)."\n");
            #}
        }
        $expression_fh->close;    

        $RNAconfig.=<<EOS;

#RNA EXPRESSION DATA 
<plot>
# The type sets the format of the track.
type = histogram
file = $output_directory/data/genes.fpkm.expsort.top1percent.tsv
min=0
max=15 #Cap, otherwise outliers hide everything

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the ideogram.
r1   = 0.80r
r0   = 0.70r

# Histograms can have both a fill and outline. The default outline is 1px thick black.
fill_color = red

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.
thickness = 0p

# Do not join histogram bins that do not abut.
extend_bin = no
orientation = out

#Draw background to highlight track
<backgrounds>
show  = yes
<background>
color = vvlgrey
</background>
</backgrounds>

#Draw axes lines
<axes>
show = yes
thickness = 1
color     = lgrey
<axis>
spacing   = 0.1r
</axis>
</axes>

</plot>
    


EOS
    }
    }else{
        $self->status_message("There was no rna data for this build. This track will be empty");
    }
    
   return $RNAconfig;     
}

sub execute {
    my $self = shift;
    
    #TODO Determine what data types should be required to run this tool
    #TODO code for wgs+exome and wgs+rna only
    
    
    # grab params from $self
    my $build = $self->build;
    my $wgs_build = $build->wgs_build;
    my $exo_build = $build->exome_build;
    
    unless($wgs_build||$exo_build){
        $self->status_message("There is no WGS or exome build for this Clinseq build. This tool cannot be run");
        return;
    }
    

    my $output_directory = $self->output_directory;

    $self->status_message("Running on build " . $build->__display_name__); 
    $self->status_message("Output directory is " . $self->output_directory);


    # initialize directories
    unless (-d $output_directory) {
        # this module has wrappers which do logging, throw exceptions, around regular tasks
        Genome::Sys->create_directory($output_directory);
    }
    unless (-d "$output_directory/data"){
        Genome::Sys->create_directory("$output_directory/data");
    }

    unless (-d "$output_directory/raw"){
        Genome::Sys->create_directory("$output_directory/raw");
    }


    my $dataDir = $build->data_directory . "/" . $build->common_name;

 

    # declare variable to contain text for circos.conf file output
    my $config;

    # Gene hash is a union of all of the different genes in each file used for the gene annotations on the plot.
    # Two gene hashes are created and two files are produced that could be used for gene labels on the circos plot
    #    1. genes-noAmpDel is a hash with all of genes except for amplifications and deletions
    #    2. genes-AmpDel is a wash with the amplifications and deletions (This file could enrich lables to only be in the AmpDel regions)
    my %genes_noAmpDel;
    my %genes_AmpDel;


    # determine if data is exome or WGS and access correct folders for input files
    
    if($wgs_build){

        ###Candidate Fusions
        Genome::Sys->copy_file("$dataDir/sv/CandidateSvCodingFusions.tsv", "$output_directory/raw/CandidateSvCodingFusions.tsv");
        ### The following if statement is necessary to integrate MakeCircosPlot into the ClinSeq pipeline.
        if($self->candidate_fusion_infile){
	        system("rm -f $output_directory/raw/CandidateSvCodingFusions.tsv");
	        Genome::Sys->copy_file($self->candidate_fusion_infile , "$output_directory/raw/CandidateSvCodingFusions.tsv");
	    }
        ###Fusions
        if(my $tumor_rnaseq_build = $build->tumor_rnaseq_build && -e $dataDir."/rnaseq/tumor/fusions/filtered_chimeras.bedpe"){
            Genome::Sys->copy_file($dataDir."/rnaseq/tumor/fusions/filtered_chimeras.bedpe", "$output_directory/raw/filtered_chimeras.bedpe");
	    }
	    ###Deletions and Focal Amplifications
        Genome::Sys->copy_file("$dataDir/clonality/cnaseq.cnvhmm", "$output_directory/raw/cnaseq.cnvhmm");
        if($self->cnv_hmm_file){
	        system("rm -f $output_directory/raw/cnaseq.cnvhmm");
	        Genome::Sys->copy_file($self->cnv_hmm_file , "$output_directory/raw/cnaseq.cnvhmm");
	    }
        ###Deletions and focal amplifications gene files
        Genome::Sys->copy_file("$dataDir/cnv/wgs_cnv/cnview/cnv.All_genes.ampdel.tsv" , "$output_directory/raw/cnv.All_genes.ampdel.tsv");
        if($self->gene_ampdel_file){
            system("rm -f $output_directory/raw/cnv.All_genes.ampdel.tsv");
            Genome::Sys->copy_file($self->gene_ampdel_file , "$output_directory/raw/cnv.All_genes.ampdel.tsv");
        }
       
        #TODO tidy this section, superfluous since exome or WGS has already been defined

        ### Tier1 SNVs and INDELs
        #decides which somatic variation model to use
        my $snv_data_dir;
        my $indel_data_dir;
        if($wgs_build && $exo_build){
            $snv_data_dir="$dataDir/snv/wgs_exome";
            $indel_data_dir="$dataDir/indel/wgs_exome";
        }elsif($exo_build){
            $snv_data_dir="$dataDir/snv/exome";
            $indel_data_dir="$dataDir/indel/exome";
        }else{
            $snv_data_dir="$dataDir/snv/wgs";
            $indel_data_dir="$dataDir/indel/wgs";
        }
        unless ($self->annotated_variants_tsv) {
            Genome::Sys->copy_file("$snv_data_dir/snvs.hq.tier1.v1.annotated.compact.tsv", "$output_directory/raw/snvs.hq.tier1.v1.annotated.compact.tsv");
            Genome::Sys->copy_file("$indel_data_dir/indels.hq.tier1.v1.annotated.compact.tsv", "$output_directory/raw/indels.hq.tier1.v1.annotated.compact.tsv");
        }
 
        $config =<<EOS;
# Chromosome name, size and color definition
karyotype = data/karyotype/karyotype.human.txt

# The <ideogram> block defines the position, size, labels and other
# properties of the segments on which data are drawn. These segments
# are usually chromosomes, but can be any integer axis.

<ideogram>

<spacing>
# Spacing between ideograms. Suffix "r" denotes a relative value. It
# is relative to circle circumference (e.g. space is 0.5% of
# circumference).
default = 0.005r
</spacing>

# Ideogram position, thickness and fill. 
#
# Radial position within the image of the ideograms. This value is
# usually relative ("r" suffix).
radius           = 0.80r

# Thickness of ideograms, which can be absolute (e.g. pixels, "p"
# suffix) or relative ("r" suffix). When relative, it is a fraction of
# image radius.
thickness        = 20p

# Ideograms can be drawn as filled, outlined, or both. When filled,
# the color will be taken from the last field in the karyotype file,
# or set by chromosomes_colors. Color names are discussed in
# http://www.circos.ca/documentation/tutorials/configuration/configuration_files
# When stroke_thickness=0p or if the parameter is missing, the ideogram is
# has no outline and the value of stroke_color is not used.

fill             = yes   # A
stroke_color     = dgrey # B
stroke_thickness = 2p    # B

# Definition for ideogram labels.
show_label       = yes
# see etc/fonts.conf for list of font names
label_font       = default 
label_radius     = 1.15r
label_size       = 30
label_parallel   = yes

#draw chromosome bands
show_bands = yes
fill_bands = yes
band_transparency = 4

</ideogram>

##Define some tick marks
#show_ticks          = yes
#show_tick_labels    = yes

#<ticks>
#radius           = 1r
#color            = black
#thickness        = 2p

## the tick label is derived by multiplying the tick position
## by 'multiplier' and casting it in 'format':
## sprintf(format,position*multiplier)
#multiplier       = 1e-6

## \%d   - integer
## \%f   - float
## \%.1f - float with one decimal
## \%.2f - float with two decimals
## for other formats, see http://perldoc.perl.org/functions/sprintf.html
#format           = \%d

#<tick>
#spacing        = 5u
#size           = 10p
#</tick>

#<tick>
#spacing        = 25u
#size           = 15p
#show_label     = yes
#label_size     = 20p
#label_offset   = 10p
#format         = \%d
#</tick>

#</ticks>
EOS

        ###Candidate Fusions
        my $candidate_fusions = Genome::Sys->read_file("$output_directory/raw/CandidateSvCodingFusions.tsv");
        my $candidate_fusions_fh = Genome::Sys->open_file_for_writing("$output_directory/data/CandidateSvCodingFusions.txt");
        while ($candidate_fusions =~ /(\S+)\s+(\S+)\s+chr(\S+):(\d+)-(\d+)\s+chr(\S+):(\d+)-(\d+)/g) {
            $genes_noAmpDel{$1}="hs$3\t$4\t$5";
            $genes_noAmpDel{$2}="hs$6\t$7\t$8";
            $genes_AmpDel{$1}="hs$3\t$4\t$5";
            $genes_AmpDel{$2}="hs$6\t$7\t$8";
            print $candidate_fusions_fh ("hs$3 $4 $5 hs$6 $7 $8\n");
        }
        $candidate_fusions_fh->close;
        $config .=<<EOS;
#FUSION DATA
<links>
<link>
file          = $output_directory/data/CandidateSvCodingFusions.txt
radius        = 0.50r
bezier_radius = 0r
color         = black_a1
thickness     = 2
</link>

EOS

        ###Fusions
        if(my $tumor_rnaseq_build = $build->tumor_rnaseq_build && -e $dataDir."/rnaseq/tumor/fusions/filtered_chimeras.bedpe"){
            my $fusions = Genome::Sys->read_file("$output_directory/raw/filtered_chimeras.bedpe");
            my $fusions_fh = Genome::Sys->open_file_for_writing("$output_directory/data/filtered_chimeras.bedpe");
            while ($fusions =~ /(\d+|X|Y)\s(\d+)\s(\d+)\s(\d+|X|Y)\s(\d+)\s(\d+)\s+\w+\s\d+\s[+|-]\s[+|-]\s(\S+):(\S+)/g) {
#           		print "gene 1 : $7 hs$1 $2 $3";
#           		print "gene 2 : $8 hs$4 $5 $6";
#                $genes_noAmpDel{$7}="hs$1 $2 $3";
#                $genes_noAmpDel{$8}="hs$4 $5 $6";
#                $genes_AmpDel{$7}="hs$1 $2 $3";
#                $genes_AmpDel{$8}="hs$4 $5 $6";
                print $fusions_fh ("hs$1 $2 $3 hs$4 $5 $6\n");
            }
            $fusions_fh->close;
            $config .=<<EOS;
#Fusions RNAseq support
<link>
file          = $output_directory/data/filtered_chimeras.bedpe
radius          = 0.50r
bezier_radius = 0r
color         = red_a1
thickness     = 2
</link>        
EOS
       }

    $config.=<<EOS;
</links>

EOS

    $config.=<<EOS;

<plots>
EOS

        ###Deletions and Focal Amplifications
        my $deletions_and_focal_amps = Genome::Sys->read_file("$output_directory/raw/cnaseq.cnvhmm");
        my $deletions_fh = Genome::Sys->open_file_for_writing("$output_directory/data/deletions.txt");
        my $focal_amps_fh = Genome::Sys->open_file_for_writing("$output_directory/data/focalAmps.txt");
        while ($deletions_and_focal_amps =~ /(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\n/g) {
        
            if($6 eq "Loss"){
                my $loss;
                if($4<$5){
                    $loss = $4-$5;
                }else{
                    $loss = $5-$4;
                }
                #this is the max deletion threshold
                if($loss>2){
                    $loss=2;
                }
                print $deletions_fh ("hs$1 $2 $3 $loss\n");
            }else{
                my $gain;
                if($4>$5){
                    $gain = $4-$5;
                }else{
                    $gain = $5-$4;
                }
                #this is the max amplification threshold
                if($gain>6){
                    $gain=6;
                }
                print $focal_amps_fh ("hs$1 $2 $3 $gain\n");
            }
        
        }
        $deletions_fh->close;
        $focal_amps_fh->close;
    
        #the gene names for deletions and focal amplifications is in another file so it needs to be read in here
        my $ampdel_genes = Genome::Sys->read_file("$output_directory/raw/cnv.All_genes.ampdel.tsv");
        while ($ampdel_genes =~ /ENS\w+\s\w+\s(\S+)\s\S+\s(\w+)\s(\d+)\s(\d+)/g) {
            $genes_AmpDel{$1}="hs$2\t$3\t$4";
        }
         


    ## if exome build data is available, collects input files from appropriate folders 

    }else{

        ###Candidate Fusions
        # not availabe for exome only data
        
        ###Fusions
        if(my $tumor_rnaseq_build = $build->tumor_rnaseq_build && -e $dataDir."/rnaseq/tumor/fusions/filtered_chimeras.bedpe"){
            Genome::Sys->copy_file($dataDir."/rnaseq/tumor/fusions/filtered_chimeras.bedpe", "$output_directory/raw/filtered_chimeras.bedpe");
        }
        ###Deletions and Focal Amplifications
        Genome::Sys->copy_file("$dataDir/cnv/exome_cnv/cnmops.cnvhmm", "$output_directory/raw/cnmops.cnvhmm");
        if($self->cnv_hmm_file){
            system("rm -f $output_directory/raw/cnmops.cnvhmm");
            Genome::Sys->copy_file($self->cnv_hmm_file , "$output_directory/raw/cnmops.cnvhmm");
        }
        ###Deletions and focal amplifications gene files
        Genome::Sys->copy_file("$dataDir/cnv/exome_cnv/CNView_All/CNView_All_genes.ampdel.tsv" , "$output_directory/raw/CNView_All_genes.ampdel.tsv");
        if($self->gene_ampdel_file){
            system("rm -f $output_directory/raw/cnv.All_genes.ampdel.tsv");
            Genome::Sys->copy_file($self->gene_ampdel_file , "$output_directory/raw/cnv.All_genes.ampdel.tsv");
        }
        
        

        #TODO tidy this section, superfluous since exome or WGS has already been defined

        ### Tier1 SNVs and INDELs
        #decides which somatic variation model to use
        my $snv_data_dir;
        my $indel_data_dir;
        if($wgs_build && $exo_build){
            $snv_data_dir="$dataDir/snv/wgs_exome";
            $indel_data_dir="$dataDir/indel/wgs_exome";
        }elsif($exo_build){
            $snv_data_dir="$dataDir/snv/exome";
            $indel_data_dir="$dataDir/indel/exome";
        }else{
            $snv_data_dir="$dataDir/snv/wgs";
            $indel_data_dir="$dataDir/indel/wgs";
        }
        unless ($self->annotated_variants_tsv) {
            Genome::Sys->copy_file("$snv_data_dir/snvs.hq.tier1.v1.annotated.compact.tsv", "$output_directory/raw/snvs.hq.tier1.v1.annotated.compact.tsv");
            Genome::Sys->copy_file("$indel_data_dir/indels.hq.tier1.v1.annotated.compact.tsv", "$output_directory/raw/indels.hq.tier1.v1.annotated.compact.tsv");
        }
        

        
    #TODO if user enters specific files to run then overide the standard files retrived from the build here
        
        $config =<<EOS;
# Chromosome name, size and color definition
karyotype = data/karyotype/karyotype.human.txt

# The <ideogram> block defines the position, size, labels and other
# properties of the segments on which data are drawn. These segments
# are usually chromosomes, but can be any integer axis.

<ideogram>

<spacing>
# Spacing between ideograms. Suffix "r" denotes a relative value. It
# is relative to circle circumference (e.g. space is 0.5% of
# circumference).
default = 0.005r
</spacing>

# Ideogram position, thickness and fill. 
#
# Radial position within the image of the ideograms. This value is
# usually relative ("r" suffix).
radius           = 0.80r

# Thickness of ideograms, which can be absolute (e.g. pixels, "p"
# suffix) or relative ("r" suffix). When relative, it is a fraction of
# image radius.
thickness        = 20p

# Ideograms can be drawn as filled, outlined, or both. When filled,
# the color will be taken from the last field in the karyotype file,
# or set by chromosomes_colors. Color names are discussed in
# http://www.circos.ca/documentation/tutorials/configuration/configuration_files
# When stroke_thickness=0p or if the parameter is missing, the ideogram is
# has no outline and the value of stroke_color is not used.

fill             = yes   # A
stroke_color     = dgrey # B
stroke_thickness = 2p    # B

# Definition for ideogram labels.
show_label       = yes
# see etc/fonts.conf for list of font names
label_font       = default 
label_radius     = 1.15r
label_size       = 30
label_parallel   = yes

#draw chromosome bands
show_bands = yes
fill_bands = yes
band_transparency = 4

</ideogram>

##Define some tick marks
#show_ticks          = yes
#show_tick_labels    = yes

#<ticks>
#radius           = 1r
#color            = black
#thickness        = 2p

## the tick label is derived by multiplying the tick position
## by 'multiplier' and casting it in 'format':
## sprintf(format,position*multiplier)
#multiplier       = 1e-6

## \%d   - integer
## \%f   - float
## \%.1f - float with one decimal
## \%.2f - float with two decimals
## for other formats, see http://perldoc.perl.org/functions/sprintf.html
#format           = \%d

#<tick>
#spacing        = 5u
#size           = 10p
#</tick>

#<tick>
#spacing        = 25u
#size           = 15p
#show_label     = yes
#label_size     = 20p
#label_offset   = 10p
#format         = \%d
#</tick>

#</ticks>
EOS


# WGS Fusions not available for exome so this will be skipped
=cut
        #TODO I am assuming that this comes from the WGS data. I am not sure on this, however if it does not the <links> tag needs to be brought out of the conf section here
        ###Candidate Fusions
        my $candidate_fusions = Genome::Sys->read_file("$output_directory/raw/CandidateSvCodingFusions.tsv");
        my $candidate_fusions_fh = Genome::Sys->open_file_for_writing("$output_directory/data/CandidateSvCodingFusions.txt");
        while ($candidate_fusions =~ /(\S+)\s+(\S+)\s+chr(\S+):(\d+)-(\d+)\s+chr(\S+):(\d+)-(\d+)/g) {
            $genes_noAmpDel{$1}="hs$3\t$4\t$5";
            $genes_noAmpDel{$2}="hs$6\t$7\t$8";
            $genes_AmpDel{$1}="hs$3\t$4\t$5";
            $genes_AmpDel{$2}="hs$6\t$7\t$8";
            print $candidate_fusions_fh ("hs$3 $4 $5 hs$6 $7 $8\n");
        }
        $candidate_fusions_fh->close;
        $config .=<<EOS;
#FUSION DATA
<links>
<link>
file          = $output_directory/data/CandidateSvCodingFusions.txt
radius        = 0.50r
bezier_radius = 0r
color         = black_a1
thickness     = 2
</link>

EOS
=cut

        ###Fusions
        if(my $tumor_rnaseq_build = $build->tumor_rnaseq_build && -e $dataDir."/rnaseq/tumor/fusions/filtered_chimeras.bedpe"){
            my $fusions = Genome::Sys->read_file("$output_directory/raw/filtered_chimeras.bedpe");
            my $fusions_fh = Genome::Sys->open_file_for_writing("$output_directory/data/filtered_chimeras.bedpe");
            while ($fusions =~ /(\d+|X|Y)\s(\d+)\s(\d+)\s(\d+|X|Y)\s(\d+)\s(\d+)\s+\w+\s\d+\s[+|-]\s[+|-]\s(\S+):(\S+)/g) {
                print $fusions_fh ("hs$1 $2 $3 hs$4 $5 $6\n");
            }
            $fusions_fh->close;
            $config .=<<EOS;
#Fusions RNAseq support
<links>
<link>
file          = $output_directory/data/filtered_chimeras.bedpe
radius          = 0.50r
bezier_radius = 0r
color         = red_a1
thickness     = 2
</link>        
EOS
       }

    $config.=<<EOS;
</links>

EOS

    $config.=<<EOS;

<plots>
EOS


        ###Deletions and Focal Amplifications
        my $deletions_and_focal_amps = Genome::Sys->read_file("$output_directory/raw/cnmops.cnvhmm");
        my $deletions_fh = Genome::Sys->open_file_for_writing("$output_directory/data/deletions.txt");
        my $focal_amps_fh = Genome::Sys->open_file_for_writing("$output_directory/data/focalAmps.txt");
        while ($deletions_and_focal_amps =~ /(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+\S+\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+(\S+)\n/g) {
        
            if($6 eq "Loss"){
                my $loss;
                if($4<$5){
                    $loss = $4-$5;
                }else{
                    $loss = $5-$4;
                }
                #this is the max deletion threshold
                if($loss>2){
                    $loss=2;
                }
                print $deletions_fh ("hs$1 $2 $3 $loss\n");
            }else{
                my $gain;
                if($4>$5){
                    $gain = $4-$5;
                }else{
                    $gain = $5-$4;
                }
                #this is the max amplification threshold
                if($gain>6){
                    $gain=6;
                }
                print $focal_amps_fh ("hs$1 $2 $3 $gain\n");
            }
        
        }
        $deletions_fh->close;
        $focal_amps_fh->close;  
    
        #the gene names for deletions and focal amplifications is in another file so it needs to be read in here
        my $ampdel_genes = Genome::Sys->read_file("$output_directory/raw/CNView_All_genes.ampdel.tsv");
        while ($ampdel_genes =~ /ENS\w+\s\w+\s(\S+)\s\S+\s(\w+)\s(\d+)\s(\d+)/g) {
            $genes_AmpDel{$1}="hs$2\t$3\t$4";
        }
        
		
    }


    ## following occurs with WGS or exome data input

    $config .=<<EOS;

#DELETIONS DATA
<plot>
# The type sets the format of the track.
type = histogram
file = $output_directory/data/deletions.txt
min=-1
max=2

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the ideogram.
r1   = 0.65r
r0   = 0.55r

# Histograms can have both a fill and outline. The default outline is 1px thick black. 
fill_color = blue 

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.
thickness = 0p

# Do not join histogram bins that do not abut.
extend_bin = no
orientation = out

#Draw background to highlight track
<backgrounds>
show  = yes
<background>
color = vvlgrey
</background>
</backgrounds>

#Draw axes lines
<axes>
show = yes
thickness = 1
color     = lgrey
<axis>
spacing   = 0.1666666666r
</axis>
<axis>
position   = 0.333333333r
color      = black
</axis>
</axes>
</plot>        

#FOCAL AMPLIFICATIONS DATA
<plot>
# The type sets the format of the track.
type = histogram
file = $output_directory/data/focalAmps.txt
min=-1
max=2

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the ideogram.
r1   = 0.65r
r0   = 0.55r

# Histograms can have both a fill and outline. The default outline is 1px thick black. 
fill_color = red

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.
thickness = 0p

# Do not join histogram bins that do not abut.
extend_bin = no
orientation = out

</plot>
EOS
	#DIFFERENTIAL EXPRESSION
	$config .= $self->rnaExpression($build, $output_directory, $dataDir, \%genes_AmpDel, \%genes_noAmpDel, $config);
        
    ### Tier1 SNVs and INDELs

    my $snv_fh = Genome::Sys->open_file_for_writing("$output_directory/data/snvs.hq.tier1.v1.annotated.compact.tsv");
    my $indel_fh = Genome::Sys->open_file_for_writing("$output_directory/data/indels.hq.tier1.v1.annotated.compact.tsv");
    
    #SNV
    unless ($self->annotated_variants_tsv ) {
        my $snv_file = Genome::Sys->read_file("$output_directory/raw/snvs.hq.tier1.v1.annotated.compact.tsv");
        while ($snv_file =~ /(\S+):(\d+)-(\d+)\s+\S+\s+(\S+)\s+\w+\s+(\S+)\s+(\S+)\s+(\S+).*/g) {
            $genes_noAmpDel{$4}="hs$1\t$2\t$3";
            $genes_AmpDel{$4}="hs$1\t$2\t$3";
            print $snv_fh "hs$1 $2 $3 0.5 fill_color=black\n";
            #        switch($5){
            #            case "nonsense"            {print $snv_fh "fill_color=goldenrod\n"}
            #            case "missense"            {print $snv_fh "fill_color=blue\n"}
            #            case "silent"            {print $snv_fh "fill_color=green\n"}
            #            case "splice_site"        {print $snv_fh "fill_color=black\n"}
            #            case "rna"                {print $snv_fh "fill_color=purple\n"}
            #            case "nonstop"            {print $snv_fh "fill_color=black\n"}
            #        }
        }

    #Indel
        my $indel_file = Genome::Sys->read_file("$output_directory/raw/indels.hq.tier1.v1.annotated.compact.tsv");
        while ($indel_file =~ /(\S+):(\d+)-(\d+)\s+\S+\s+(\S+)\s+\w+\s+(\S+)\s+(\S+)\s+(\S+).*/g) {
            $genes_noAmpDel{$4}="hs$1\t$2\t$3";
            $genes_AmpDel{$4}="hs$1\t$2\t$3";
            my $color;
            if ($6 eq "-"){$color="red";}
            else {$color="blue";}
            print $indel_fh "hs$1 $2 $3 0.5 fill_color=$color\n";
        }
    } else {
        my $var_file = Genome::Utility::IO::SeparatedValueReader->create(
            input => $self->annotated_variants_tsv,
        );
        while (my $data = $var_file->next) {
            $genes_noAmpDel{$data->{'default_gene_name'}} = 'hs'. $data->{'chromosome'} ."\t". $data->{'start'} ."\t". $data->{'stop'};
            $genes_AmpDel{$data->{'default_gene_name'}} = 'hs'. $data->{'chromosome'} ."\t". $data->{'start'} ."\t". $data->{'stop'};
            if ($data->{'type'} eq 'snv') {
                print $snv_fh 'hs'. $data->{'chromosome'} .' '. $data->{'start'} .' '. $data->{'stop'} ." 0.5 fill_color=black\n";
            } elsif ($data->{'type'} eq 'indel') {
                my $color;
                if ($data->{'reference'} eq "-") {
                    # Insertion
                    $color="red";
                }  else {
                    $color="blue";
                }
                print $indel_fh 'hs'. $data->{'chromosome'} .' '. $data->{'start'} .' '. $data->{'stop'} ." 0.5 fill_color=$color\n";;
            }
        }
    }
    $indel_fh->close;
    $snv_fh->close;

    $config .=<<EOS;
#TIER 1 SNV DATA
<plot>
# The type sets the format of the track.
type = histogram
file = $output_directory/data/snvs.hq.tier1.v1.annotated.compact.tsv
min=0
max=0.5

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the ideogram.
r1   = 0.95r
r0   = 0.85r

# Histograms can have both a fill and outline. The default outline is 1px thick black.
#fill_color = blue  #color specified as option in data file

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.
thickness = 0p

# Do not join histogram bins that do not abut.
extend_bin = no
orientation = out

#Draw background to highlight track
<backgrounds>
show  = yes
<background>
color = vvlgrey
</background>
</backgrounds>
</plot>



#TIER 1 INDEL DATA
<plot>
# The type sets the format of the track.
type = histogram
file = $output_directory/data/indels.hq.tier1.v1.annotated.compact.tsv
min=0
max=0.5

# The track is confined within r0/r1 radius limits. When using the
# relative "r" suffix, the values are relative to the position of the ideogram.
r1   = 0.95r
r0   = 0.85r

# Histograms can have both a fill and outline. The default outline is 1px thick black.
fill_color = black

# To turn off default outline, set the outline thickness to zero. If
# you want to permanently disable this default, edit
# etc/tracks/histogram.conf in the Circos distribution.
thickness = 0p

# Do not join histogram bins that do not abut.
extend_bin = no
orientation = out
</plot>        
EOS
    ###Gene List
    my $gene_fh = Genome::Sys->open_file_for_writing("$output_directory/raw/genes_noAmpDel.txt");
    print $gene_fh "chr\tstart\tend\tgene\n";
    foreach my $gene(keys %genes_noAmpDel){
        print $gene_fh "$genes_noAmpDel{$gene}\t$gene\n";
    }
    $gene_fh->close;
	if(!$build->cancer_annotation_db){
	    Genome::Sys->status_message("ERROR: No cancer_annotation_db");
	    return 1;
	}
    my $cancer_annotation_db = $build->cancer_annotation_db;
    #my $cancer_annotation_db = Genome::Db->get(id => 'tgi/cancer-annotation/human/build37-20130711.1');
    my $annotate_genes_cmd1 = Genome::Model::ClinSeq::Command::AnnotateGenesByCategory->create(
        infile => "$output_directory/raw/genes_noAmpDel.txt",
        cancer_annotation_db => $cancer_annotation_db,
        gene_name_columns => ['gene'],
    );
    $annotate_genes_cmd1->execute() or die;
    
    open (ANNOT_FILE, "<$output_directory/raw/genes_noAmpDel.catanno.txt") or die "\nCouldn't open genes_noAmpDel.catanno.txt file.";
    my $head = <ANNOT_FILE>;
    my @columns = split('\t',$head);
    my $last_column = scalar @columns;
    close (ANNOT_FILE);

    my $sort_cmd1 = "sort -rnk $last_column $output_directory/raw/genes_noAmpDel.catanno.txt|head -100 |cut -d \"\t\" -f 1-4  > $output_directory/data/genes_noAmpDel.catanno.sorted.txt";
    Genome::Sys->shellcmd(cmd => $sort_cmd1, set_pipefail => 1,);



    my $geneAmpDel_fh = Genome::Sys->open_file_for_writing("$output_directory/raw/genes_AmpDel.txt");
    print $geneAmpDel_fh "chr\tstart\tend\tgene\n";
    foreach my $gene(keys %genes_AmpDel){
        print $geneAmpDel_fh "$genes_AmpDel{$gene}\t$gene\n";
    }
    $geneAmpDel_fh->close;

    my $annotate_genes_cmd2 = Genome::Model::ClinSeq::Command::AnnotateGenesByCategory->create(
        infile => "$output_directory/raw/genes_AmpDel.txt",
        cancer_annotation_db => $cancer_annotation_db,
        gene_name_columns => ['gene'],
    );
    $annotate_genes_cmd2->execute() or die;
    my $sort_cmd2 = "sort -rnk $last_column $output_directory/raw/genes_AmpDel.catanno.txt|head -100 |cut -d \"\t\" -f 1-4  > $output_directory/data/genes_AmpDel.catanno.sorted.txt";
    Genome::Sys->shellcmd(cmd => $sort_cmd2, set_pipefail => 1,);

    $config .=<<EOS;
#GENE LABELS
<plot>
type  = text
file  = $output_directory/data/genes_noAmpDel.catanno.sorted.txt

# Like with other tracks, text is limited to a radial range by setting
# r0 and r1.
#
# Individual labels can be repositioned automatically with in a
# position window to fit more labels, without overlap. This is an
# advanced feature - see the 2D Track text tutorials.
r1    = 1.1r
r0    = 1.0r

# For a list of fonts, see etc/fonts.conf in the Circos distribution.
#label_font = light
label_font = condensed
label_size = 20p
#label_rotate = yes

# padding  - text margin in angular direction
# rpadding - text margin in radial direction
rpadding   = 5p

# Short lines can be placed before the label to connect them to the
# label's position. This is most useful when the labels are
# rearranged.
show_links     = yes
link_dims      = 0p,2p,5p,2p,2p
link_thickness = 2p
link_color     = black

#Turn on "snuggling" to allow more labels to fit
label_snuggle         = yes
max_snuggle_distance  = 1r
snuggle_tolerance     = 0.25r
snuggle_sampling      = 2
snuggle_link_overlap_test = yes
snuggle_link_overlap_tolerance = 2p
snuggle_refine        = yes

</plot>
    
EOS



=cut    
    # example accessing data:
    # get an input from the build
    my $wgs_build = $build->wgs_build;

    TODO add error checking
    # all inputs are optional, so test for it being there before using it
    if ($wgs_build) {
        $self->status_message("This clinseq build has WGS somatic data, build: " . $wgs_build->__display_name__);
    
        # when you have a build, you can get the path for a file in the data directory through an API call
        my $wgs_somatic_snvs_tier1_path = $wgs_build->data_set_path('effects/snvs.hq.novel.tier1','2','bed');
        $self->status_message("WGS somatic variants are at $wgs_somatic_snvs_tier1_path");
    
        # the above is the same as the following, but the following circumvents the API
        #my $wgs_somatic_snvs_tier1_path = $wgs_build->data_directory . '/effects/snvs.hq.novel.tier1.v2.bed'; 
    }
    else {
        $self->status_message("No WGS data on this clinseq model.");
    }

    # get an input from an input
    my $wgs_tumor_refalign = $wgs_build->tumor_build;
=cut




    $config.=<<EOS;

</plots>

<colors>
goldenrod = 218,165,32
teal = 180,100,25
chocolate = 210,105,30
</colors>

################################################################
# The remaining content is standard and required. It is imported from
# default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files, 
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>                
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>> 

# Debugging, I/O an dother system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>> 
EOS
    ## write the config file to point to the new data files above
    my $out_fh = Genome::Sys->open_file_for_writing("$output_directory/circos.conf");
    print $out_fh $config;
    $out_fh->close;

    # run circos using the correct path for the specified version
    # errors will throw an exception

    my $path_to_circos_executable = Genome::Sys->sw_path("circos",$self->use_version);
    $self->status_message("using path $path_to_circos_executable");


   Genome::Sys->shellcmd(cmd => "cd $output_directory; $path_to_circos_executable -conf ./circos.conf");


    return 1;
}

1;

