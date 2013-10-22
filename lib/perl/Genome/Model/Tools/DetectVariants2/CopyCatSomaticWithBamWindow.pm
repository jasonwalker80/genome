package Genome::Model::Tools::DetectVariants2::CopyCatSomaticWithBamWindow;

use warnings;
use strict;

use Cwd;
use Genome;
use Workflow::Simple;

my $DEFAULT_VERSION = '0.1';

class Genome::Model::Tools::DetectVariants2::CopyCatSomaticWithBamWindow{
    is => ['Genome::Model::Tools::DetectVariants2::Detector'],
    doc => "Produces somatic copy-number calls from paired samples",
    has_param => [
        lsf_resource => {
            default_value => 'rusage[mem=4000] select[type==LINUX64 && maxtmp>10000] span[hosts=1]',
        },
    ],
};


sub _detect_variants {
    my $self = shift;

    ##parse input params string - expected format
    #--bamwindow-version 0.4 --bamwindow-params [-w 10000 -r -l -s]
    my $params = $self->params;
    my $bamwindow_version;
    my $bamwindow_params;

    if ($params =~ m/--bamwindow-version/) {
        $params =~ m/--bamwindow-version\s*(\d+\.?\d?)\s*/;
        $bamwindow_version = $1;
        $params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
        $bamwindow_params  = $1;
        $bamwindow_params =~ s/--bamwindow-params\s*//;
    }

    my $copycat_params;
    #--copycat-params [--per-read-length --per-library]
    if ($params =~ m/--copycat-params/) {
        $params =~ m/--bamwindow-params\s*\{([^\}]+)\}/;
        $bamwindow_params = $1;
        $params =~ m/--copycat-params\s*\{([^\}]+)\}/;
        $copycat_params = $1;
    }
    my $per_library = 0;
    my $per_read_length = 0;
    if($copycat_params =~ m/--per-library/) {
        $per_library = 1;
    }
    if($copycat_params =~ m/--per-read-length/) {
        $per_read_length = 1;
    }



    my %input;

    # Define a workflow from the static XML at the bottom of this module
    my $workflow = Workflow::Operation->create_from_xml(\*DATA);

    # Validate the workflow
    my @errors = $workflow->validate;
    if (@errors) {
        $self->error_message(@errors);
        die "Errors validating workflow\n";
    }

    my $genome_build = $self->_get_genome_build; #TODO: This is literally the dumbest thing on earth.  Fix it

    # Collect and set input parameters
    #for bam-window
    $input{bamwindow_version} = $bamwindow_version;
    $input{bamwindow_params} = $bamwindow_params;
    $input{tumor_bam_file} = $self->aligned_reads_input;
    $input{normal_bam_file} = $self->control_aligned_reads_input;
    $input{tumor_window_file} = join('/', $self->_temp_staging_directory, 'tumor_window_file');
    $input{normal_window_file} = join('/', $self->_temp_staging_directory, 'normal_window_file');

    # my $normal_window_dir = join('/', $self->_temp_staging_directory, 'normal_window_dir');
    # mkdir $normal_window_dir;
    # $input{normal_window_dir} = $normal_window_dir;
    # my $tumor_window_dir = join('/', $self->_temp_staging_directory, 'tumor_window_dir');
    # mkdir $tumor_window_dir;
    # $input{tumor_window_dir} = $tumor_window_dir; 

    #for copycat
    $input{per_read_length} = $per_read_length;
    $input{per_library} = $per_library;
    $input{tumor_samtools_file} = $self->get_samtools_results($self->aligned_reads_input);
    $input{normal_samtools_file} = $self->get_samtools_results($self->control_aligned_reads_input);
    $input{copycat_output_directory} = $self->_temp_staging_directory;
    $input{annotation_directory} = '/gscmnt/gc6122/info/medseq/annotations/copyCat'; #TODO: don't do this, this is really bad
    $input{genome_build} = $genome_build;


    # $self->_dump_workflow($workflow);
    print "DEBUG: INPUTS:", Data::Dumper->Dumper(%input), "\n";

    # my $log_dir = $self->output_directory;
    my $log_dir = '/tmp/awesome';
    # if(Workflow::Model->parent_workflow_log_dir) {
        # $log_dir = Workflow::Model->parent_workflow_log_dir;
    # }
    $workflow->log_dir($log_dir);

    # Launch workflow
    $self->status_message("Launching workflow now.");
    my $result = Workflow::Simple::run_workflow_lsf($workflow, %input);

    # Collect and analyze results
    unless($result){
        if (@Workflow::Simple::ERROR){
            print Data::Dumper->Dumper(@Workflow::Simple::ERROR), "\n";
        }
        die $self->error_message("Workflow did not return correctly");
    }
    $self->_workflow_result($result);

    return 1;
}



sub has_version {
    return 1; #FIXME implement this when this module is filled out
}

sub _sort_detector_output {
    return 0;
}

##todo comment
sub get_samtools_results{
    my $self = shift;
    my $bam_path = shift;
    my $return_snp_file;
    #TODO: This will never work in testing for obvious reasons
    my @results = Genome::Model::Tools::DetectVariants2::Result::DetectionBase->get(
        detector_name => "Genome::Model::Tools::DetectVariants2::Samtools",
        aligned_reads =>$bam_path );
    if(@results) {
        return $results[0]->path("snvs.hq");
    } else {
        $self->status_message("Could not find any DV2::Samtools result object for $bam_path");
        #alternative lookup - maybe later?
    }
    return "";
}

sub _get_genome_build{
    my $self = shift;
    my $reference_build = Genome::Model::Build->get($self->reference_build_id);
    return 37 if $reference_build->build_name =~ /37/;
    return 'mm9' if $reference_build->build_name =~/mm9/;
    return 36;
}

1;

__DATA__
<?xml version='1.0' standalone='yes'?>
    <workflow name="CopyCatSomatic Detect Variants Module">


    <link fromOperation="input connector" fromProperty="normal_bam_file" toOperation="BamWindow Normal" toProperty="bam_file"/>
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Normal" toProperty="options" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Normal" toProperty="version" />
    <link fromOperation="input connector" fromProperty="normal_window_file" toOperation="BamWindow Normal" toProperty="output_file" />
    <!-- <link fromOperation="input connector" fromProperty="reference_build_id" toOperation="BamWindow Normal" toProperty="reference_build_id" /> -->


    <link fromOperation="input connector" fromProperty="tumor_bam_file" toOperation="BamWindow Tumor" toProperty="bam_file" />
    <link fromOperation="input connector" fromProperty="bamwindow_params" toOperation="BamWindow Tumor" toProperty="options" />
    <link fromOperation="input connector" fromProperty="bamwindow_version" toOperation="BamWindow Tumor" toProperty="version" />
    <link fromOperation="input connector" fromProperty="tumor_window_file" toOperation="BamWindow Tumor" toProperty="output_file" />
    <!-- <link fromOperation="input connector" fromProperty="reference_build_id" toOperation="BamWindow Tumor" toProperty="reference_build_id" /> -->


    <link fromOperation="BamWindow Normal" fromProperty="output_file" toOperation="CopyCat Somatic" toProperty="normal_window_file" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_file" toOperation="CopyCat Somatic" toProperty="tumor_window_file" />
    <link fromOperation="input connector" fromProperty="per_library" toOperation="CopyCat Somatic" toProperty="per_library" />
    <link fromOperation="input connector" fromProperty="per_read_length" toOperation="CopyCat Somatic" toProperty="per_read_length" />
    <link fromOperation="input connector" fromProperty="tumor_samtools_file" toOperation="CopyCat Somatic" toProperty="tumor_samtools_file" />
    <link fromOperation="input connector" fromProperty="normal_samtools_file" toOperation="CopyCat Somatic" toProperty="normal_samtools_file" />
    <!-- <link fromOperation="input connector" fromProperty="reference_build_id" toOperation="CopyCat Somatic" toProperty="reference_build_id" /> -->
    <!-- <link fromOperation="input connector" fromProperty="copycat_params" toOperation="CopyCat Somatic" toProperty="params" /> --> 
    <link fromOperation="input connector" fromProperty="annotation_directory" toOperation="CopyCat Somatic" toProperty="annotation_directory" />
    <link fromOperation="input connector" fromProperty="genome_build" toOperation="CopyCat Somatic" toProperty="genome_build" />
    <link fromOperation="input connector" fromProperty="copycat_output_directory" toOperation="CopyCat Somatic" toProperty="output_directory" />


    <link fromOperation="BamWindow Normal" fromProperty="output_file" toOperation="output connector" toProperty="bam_window_normal_output_file" />
    <link fromOperation="BamWindow Tumor" fromProperty="output_file" toOperation="output connector" toProperty="bam_window_tumor_output_file" />
    <link fromOperation="CopyCat Somatic" fromProperty="output_directory" toOperation="output connector" toProperty="copycat_output_directory" />

    <operation name="BamWindow Normal">
    <operationtype commandClass="Genome::Model::Tools::BamWindow" typeClass="Workflow::OperationType::Command" />
    </operation>

    <operation name="BamWindow Tumor">
    <operationtype commandClass="Genome::Model::Tools::BamWindow" typeClass="Workflow::OperationType::Command" />
    </operation>

    <operation name="CopyCat Somatic">
    <operationtype commandClass="Genome::Model::Tools::CopyCat::Somatic" typeClass="Workflow::OperationType::Command" />
    </operation>


    <operationtype typeClass="Workflow::OperationType::Model">
    <inputproperty>tumor_bam_file</inputproperty>
    <inputproperty>normal_bam_file</inputproperty>
    <inputproperty>tumor_window_file</inputproperty>
    <inputproperty>normal_window_file</inputproperty>
    <inputproperty>bamwindow_params</inputproperty>
    <inputproperty>bamwindow_version</inputproperty>
    <!-- <inputproperty>reference_build_id</inputproperty> --> 
    <inputproperty>per_read_length</inputproperty>
    <inputproperty>per_library</inputproperty>
    <inputproperty>tumor_samtools_file</inputproperty>
    <inputproperty>normal_samtools_file</inputproperty>
    <!-- <inputproperty>copycat_params</inputproperty> -->
    <inputproperty>copycat_output_directory</inputproperty>
    <inputproperty>annotation_directory</inputproperty>
    <inputproperty>genome_build</inputproperty>

    <outputproperty>bam_window_normal_output_file</outputproperty>
    <outputproperty>bam_window_tumor_output_file</outputproperty>
    <outputproperty>copycat_output_directory</outputproperty>
    </operationtype>

    </workflow>
