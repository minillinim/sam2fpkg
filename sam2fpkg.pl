#!/usr/bin/perl
###############################################################################
#
#    sam2fpkg.pl
#    
#    Converts a sam file to fpkg measurements - made for bacterial genomes
#    with no introns.
#
#    Copyright (C) 2011,2012 Michael Imelfort
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#CPAN modules

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params and print copyright
printAtStart();
my $options = checkParams();

######################################################################
# CODE HERE
######################################################################
# globals
my $global_out_file_name = $options->{'sam'}.".fpkg";
if(exists($options->{'out'}))
{
   $global_out_file_name = $options->{'out'};
}
open my $out_fh, ">", $global_out_file_name or die "**Error: could not open output file $!\n";

# contigs
my %global_contig_headers = ();             # contigs we care about
my %global_contig_lengths = ();

# counters
my $global_gene_count = 0;
my $global_non_gene_count = 0;
my $global_total_fragments_mapped = 0;
my $global_total_ORF_fragments_mapped = 0;

# string mapping
my $global_string_ID_counter = 0;           # we would prefer to store every string only once (twice at most...)
my %global_string_2_ID_map = ();
my %global_ID_2_string_map = ();

# contig to mapping possies
my %global_gene_regions = ();               # map of contig headers to start / stop positions
my %global_non_gene_regions = ();           # map of contig headers to non orf start / stop positions

# ORF maps
my %global_ORF_map = ();                    # each ORF looks like "CONTIG_ID,start_stop"
my %global_ORF_frag_counts_map = ();        # number of fragments mapped for each ORF

#non_ORF maps (non-mapped-orfs
my %global_non_ORF_map = ();
my %global_non_ORF_frag_counts_map = ();

# first open the reference file and work out which contigs we're mapping against
print "Parsing reference sequence file: ". $options->{'ref'} . "\n";
open my $ref_fh, "<", $options->{'ref'} or die "**ERROR: could not open reference file $!\n";

my $seq = "";
my $header_ID = 0;
while(<$ref_fh>)
{
    chomp $_;
    if($_ =~ /^>/)
    {
        # header line - clean if up
        $_ =~ s/^>//;
        print "\t adding contig: $_\n";
        
        if($header_ID != 0)
        {
            # all we want is the sequence length
            $global_contig_lengths{$header_ID} = length($seq);
        }
        
        # add the string to the string map and the header map
        $header_ID = addString($_); 
        $global_contig_headers{$header_ID} = 1;
        
        # make the map to hold start stops VS ORFs map for this guy
        my %tmp_map = ();
        $global_gene_regions{$header_ID} = \%tmp_map;
        
        $seq = "";
        
    }
    else
    {
        $seq .= $_;
    }
}
# last one!
if($header_ID != 0)
{
    # all we want is the sequence length
    $global_contig_lengths{$header_ID} = length($seq);
}

close $ref_fh;

# now open the gff file and work out the gene boundaries on these contigs
print "Parsing gff file: ". $options->{'gff'} . "\n";
open my $gff_fh, "<", $options->{'gff'} or die "**ERROR: could not open gff file $!\n";

# use these to work out the unmapped regions
my $last_orf_end = 1;
my $last_con_ID = 0;
while(<$gff_fh>)
{
    next if ($_ =~ /^#/);
    chomp $_;
    my @gff_fields = split /\t/, $_;
    
    my $con_ID = getID($gff_fields[0]); 
    if($last_con_ID != 0 and $con_ID != $last_con_ID and exists $global_contig_headers{$last_con_ID})
    {
        # conitg has changed!
        if($last_orf_end < $global_contig_lengths{$last_con_ID})
        {
           # addNewNonORF($last_con_ID, $last_orf_end, $global_contig_lengths{$last_con_ID});
           addNewNonORF($gff_fields[0], $last_orf_end, $global_contig_lengths{$last_con_ID});
        } 
        
        #reset this guy
        $last_orf_end = 1;
        $last_con_ID = 0;
        print "\n";
    }

    if(exists $global_contig_headers{$con_ID})
    {
        # save this guy    
        $last_con_ID = $con_ID;
        
        # we've seen this mo-fo before
        my $orf_start = int($gff_fields[3]);
        my $orf_end = int($gff_fields[4]);

        if($orf_start > $last_orf_end)
        {
            # we should chuck in a non-orf here
            # addNewNonORF($con_ID, $last_orf_end, $orf_start);
            addNewNonORF($gff_fields[0], $last_orf_end, $orf_start);
        }
        $last_orf_end = $orf_end + 1;
        
        # make a new ORF ID
        # addNewORF($con_ID, $orf_start, $orf_end);
        addNewORF($gff_fields[0], $orf_start, $orf_end);
    }
}

# do the last guy
if($last_con_ID != 0 and exists $global_contig_headers{$last_con_ID})
{
    if($last_orf_end < $global_contig_lengths{$last_con_ID})
    {
        # addNewNonORF($last_con_ID, $last_orf_end, $global_contig_lengths{$last_con_ID});
        addNewNonORF(getString($last_con_ID), $last_orf_end, $global_contig_lengths{$last_con_ID});
    }
}
 
close $gff_fh;

print "Identified: $global_gene_count genes and $global_non_gene_count unmapped regions\n";

# finally open the sam file and work out which reads hit where
print "Parsing sam file: ". $options->{'sam'} . "\n";
open my $sam_fh, "<", $options->{'sam'} or die "**ERROR: could not open sam file $!\n";
while(<$sam_fh>)
{
    next if ($_ =~ /^@/);
    chomp $_;
    my @sam_fields = split(/\t/, $_);
    my @mapping_flags = split( //, dec2bin($sam_fields[1]));
    
    # make sure it mapped
    next if($mapping_flags[4] eq "1");
    
    # get the contig ID
    # my $con_ID = getID($sam_fields[2]);
    my $con_ID = $sam_fields[2];

    if(0 != getID($con_ID))
    {
        # this is one of our guys!
        # find the orf ID
        if(exists ${$global_gene_regions{$con_ID}}{int($sam_fields[3])})
        {
            # read mapped into an orf
            my $orf_ID = ${$global_gene_regions{$con_ID}}{int($sam_fields[3])};
            
            # increment the local counter
            $global_ORF_frag_counts_map{$orf_ID}++;
            
            # increment the ORF spceific counter
            $global_total_ORF_fragments_mapped++;
        }
        elsif(exists ${$global_non_gene_regions{$con_ID}}{int($sam_fields[3])})
        {
            # read mapped elsewhere...
            my $non_orf_ID = ${$global_non_gene_regions{$con_ID}}{int($sam_fields[3])};
        }
        
        # up the total number mapped
        $global_total_fragments_mapped++;
    }
}
close $sam_fh;

# print out a heap of info
print "Total mapped: $global_total_fragments_mapped\n";
print "Mapped to ORFs: $global_total_ORF_fragments_mapped\n";
print $out_fh "#***************************\n";
print $out_fh "# ORFS\n";
print $out_fh "#***************************\n";
print $out_fh  "#ORFID\tFPKM\tLEN\tNUM\n";
foreach my $orf_ID (sort ORFSort (keys %global_ORF_map))
{
    my $orf_string = getString($orf_ID);
    my $num_mapped_to_orf = $global_ORF_frag_counts_map{$orf_ID};
    my $orf_length = $global_ORF_map{$orf_ID};
    my $fpkg = ((($num_mapped_to_orf * 1000) / $orf_length) * 1000000) / $global_total_fragments_mapped;
    print $out_fh "$orf_string\t$fpkg\t$orf_length\t$num_mapped_to_orf\n";
}
print $out_fh "#***************************\n";
print $out_fh "# NON ORFS\n";
print $out_fh "#***************************\n";
print $out_fh  "#ORFID\tFPKM\tLEN\tNUM\n";
foreach my $orf_ID (sort ORFSort (keys %global_non_ORF_map))
{
    my $orf_string = getString($orf_ID);
    my $num_mapped_to_orf = $global_non_ORF_frag_counts_map{$orf_ID};
    my $orf_length = $global_non_ORF_map{$orf_ID};
    my $fpkg = ((($num_mapped_to_orf * 1000) / $orf_length) * 1000000) / $global_total_fragments_mapped;
    print $out_fh "$orf_string\t$fpkg\t$orf_length\t$num_mapped_to_orf\n";
}

close $out_fh;

######################################################################
# CUSTOM SUBS
######################################################################
sub dec2bin { my ($dec) = @_; return sprintf "%08b", $dec; }

sub ORFSort {
    #-----
    # sort ORFs by contig and position
    #
   my @a_fields = split(/,/, getString($a)); 
   my @b_fields = split(/,/, getString($b));
   
   if(int(getID($a_fields[0])) > int(getID($b_fields[0])))
   {
       return 1;
   }  
   elsif (int(getID($a_fields[0])) == int(getID($b_fields[0])))
   {
       return int($a_fields[1]) > int($b_fields[1]);
   } 
}

sub addNewORF
{
    #-----
    # add a new ORF to the global maps
    #
    my ($con_ID, $orf_start, $orf_end) = @_;
    # get a new ID
    my $orf_str_ID = sprintf("%s,%d,%d",$con_ID, $orf_start, $orf_end);
    my $orf_ID = addString($orf_str_ID);
    # store the length of the orf
    $global_ORF_map{$orf_ID} = $orf_end - $orf_start + 1;
    $global_ORF_frag_counts_map{$orf_ID} = 0;
    # add to the lookup table so we can
    foreach my $i ($orf_start..$orf_end)
    {
        ${$global_gene_regions{$con_ID}}{$i} = $orf_ID;
    }
    # increment the total number of genes seen
    $global_gene_count++;
    print "+";
}

sub addNewNonORF
{
    #-----
    # add a new non-ORF region to the maps
    #
    my ($con_ID, $last_orf_end, $orf_start) = @_;
    # get a new ID
    my $non_orf_str_ID = sprintf("%s,%d,%d",$con_ID, $last_orf_end, $orf_start);
    my $non_orf_ID = addString($non_orf_str_ID);
    # store the length of the orf
    $global_non_ORF_map{$non_orf_ID} = $orf_start - $last_orf_end + 1;
    $global_non_ORF_frag_counts_map{$non_orf_ID} = 0;
    # add to the lookup table so we can
    foreach my $i ($last_orf_end..$orf_start)
    {
        ${$global_non_gene_regions{$con_ID}}{$i} = $non_orf_ID;
    }
    $global_non_gene_count++;
    print "-";
}

sub addString
{
    #-----
    # add a string to the stringmap
    # 
    my ($string) = @_;
    if(exists $global_string_2_ID_map{$string})
    {
        print "Warning: adding \"$string\" twice!\n";
        return $global_string_2_ID_map{$string};
    }
    $global_string_ID_counter++;
    $global_string_2_ID_map{$string} = $global_string_ID_counter;
    $global_ID_2_string_map{$global_string_ID_counter} = $string;
    return $global_string_ID_counter;
}

sub getID
{
    #-----
    # given an ID return the string (or not)
    #
    my ($string) = @_;
    if(exists $global_string_2_ID_map{$string})
    {
        return $global_string_2_ID_map{$string};
    }
    return 0;
}

sub getString
{
    #-----
    # given an ID return the string (or not)
    #
    my ($id) = @_;
    if(exists $global_ID_2_string_map{$id})
    {
        return $global_ID_2_string_map{$id};
    }
    print "Unknown ID $id\n";
    return "UNSET";
}

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    my @standard_options = ( "help|h+", "sam|s:s", "gff|g:s", "ref|r:s", "out|o:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    if(!exists $options{'sam'} ) { print "**ERROR: You need to supply a sam file\n"; exec("pod2usage $0"); }
    if(!exists $options{'gff'} ) { print "**ERROR: You need to supply a gff file\n"; exec("pod2usage $0"); }
    if(!exists $options{'ref'} ) { print "**ERROR: You need to supply a reference sequence file\n"; exec("pod2usage $0"); }
    #if(!exists $options{''} ) { print "**ERROR: \n"; exec("pod2usage $0"); }

    return \%options;
}

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) 2011,2012 Michael Imelfort
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    sam2fpkg.pl

=head1 COPYRIGHT

   copyright (C) 2011,2012 Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    Converts a sam file to fpkg measurements

=head1 SYNOPSIS

    sam2fpkg.pl -s SAMFILE -r REFFILE -g GFF3FILE [-o OUTFILE] [-help|h]

      -sam|s SAMFILE              Sam file to parse 
      -ref|r REFFILE              Reference sequence the reads were mapped to 
      -gff|g GFF3FILE             Gff3 file to parse
      -out|o OUTFILE              File to print results to (default: SAMFILE.fpkg)
      [-help -h]                  Displays basic usage information
         
=cut
