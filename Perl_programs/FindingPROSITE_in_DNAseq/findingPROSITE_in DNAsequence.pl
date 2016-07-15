#!/usr/bin/perl -w


#REMARKS: 1. In this program "print" command is used for the output of the
#         data we expect. In that case we can use it as a CGI script or
#         simply redirect the output to a file.
#         2. On the other hand "die" is used for error exits. It generates
#         messages into Error Output.

#use strict;

my %fasta_hash;      # global hash (associative array) with FASTA data
my @nsest_array;     # global array with NSEST data
my @result_array;    # global array with pattern matching results
my $pattern;         # global scalar with PROSITE pattern

# the whole program consists of few steps
# step 1. get the input parameters
my($nsest, $fasta, $prosite) = getArguments(@ARGV);

# step 2. read input data from files
readFiles($nsest, $fasta, $prosite); # could be readFiles(getArguments(@ARGV));

# step 3. check if the input data is consistent
checkDataConsistency();

# step 4. find the sequence for all exons in the get the @nsest_array
findExonsSequences();

# step 5. find all exons that match the PROSITE pattern
findExonsMatchingPattern($pattern);

# step 6. display data in the HTML form
my @result_headers = ('Name', 'Start', 'End', 'Protein');
unshift(@result_array, \@result_headers);
htmlPrintHeader("Protein search results");
htmlPrintParagraph("This table contains matches of sequence $pattern");
htmlPrintTable(@result_array);
htmlPrintFooter();


#===============================================================================
sub getArguments
#===============================================================================
{
    my(@args) = @_;
    # get command-line arguments from @ARGV array
    # if there are no arguments then print USAGE and exit
    my $USAGE = "This program requires three arguments:\n$0 <NSEST file> <FASTA file> <PROSITE file>\n\n";
    if ((@args + 0) != 3) {die $USAGE;}
    return ($args[0], $args[1], $args[2]);
}

#===============================================================================
sub readFiles
#===============================================================================
{
    my(@params) = @_;
    for (my $count = 0; $count < 3; $count++)
    {
        # does the file exist?
        unless (-e $params[$count])
        {
            die "ERROR: File ", $params[$count], " doesn\'t exist!\n";
        }
        if ($count == 0) {readNSESTFile($params[$count]);}
        elsif ($count == 1) {readFASTAFile($params[$count]);}
        elsif ($count == 2) {readPROSITEFile($params[$count]);}
    } # for
    return 0; # success
}


#===============================================================================
sub readNSESTFile
#===============================================================================
{
    my $filename = $_[0]; # first parameter
    # can we open the file?
    unless (open(NSESTFILE, $filename)) {die "ERROR: Cannot open file ", $filename, "!\n";}

    while (<NSESTFILE>)
    {
        chomp;
        if ($_ ne "") # skip empty lines
        {
            # there are actually 5 fields in the NSEST line,
            # but we add one extra for future pattern matching
            my @nsest_line = split / /, $_, 5+1;
            push @nsest_array, \@nsest_line;
        }
    } # while

    close NSESTFILE;
    return 0; # success
}

#===============================================================================
sub readFASTAFile
#===============================================================================
{
    my $filename = $_[0]; # first parameter
    # can we open the file?
    unless (open(FASTAFILE, $filename))
    {
        die "ERROR: Cannot open file ", $filename, "!\n";
    }
    my(@names, @seqs);
    my $count_names = -1;
    while (<FASTAFILE>)
    {
        chomp;
        if (/>/)
        {
            # read headers
            s/>//; # remove the starting character
            $count_names++;
            $names[$count_names] = $_;
            $seqs[$count_names] = "";
        }
        else
        {
            # read sequences in upper case
            $seqs[$count_names] .= uc($_);
        }
    } # while

    close FASTAFILE;
    #------------------------------------------
    # run a four-step check on the FASTA file
    for (my $i = 0; $i <= $count_names; $i++)
    {
        # step 1 - look for empty names
        if ($names[$i] eq "")
        {
            die "ERROR: FASTA file ", $filename, " contains empty sequence name!\n";
        }
        # step 2 - look for empty sequences
        if ($seqs[$i] eq "")
        {
            die "ERROR: FASTA file ", $filename, " contains empty sequence for ", $names[$i], "!\n";
        }
        # step 3 - look for sequences with characters other than
        if ($seqs[$i] =~ /[^ACGT]/)
        {
            die "ERROR: FASTA file ", $filename, " contains sequence with characters other than ACGT (", $names[$i], ")!\n";
        }
    }
    # step 4 - look for name duplicates
    my @sorted_names = sort {$a cmp $b} @names;
    for (my $i = 1; $i <= $count_names; $i++)
    {
        if ($sorted_names[$i - 1] eq $sorted_names[$i])
        {
            die "ERROR: FASTA file ", $filename, " contains duplicate headers ", $sorted_names[$i], "!\n";
        }
    }
    # now everything is OK, so fill the hash array
    for (my $i = 0; $i <= $count_names; $i++)
    {
        $fasta_hash{$names[$i]} = $seqs[$i];
    }
    return 0; # success
}

#===============================================================================
sub readPROSITEFile
#===============================================================================
{
    my $filename = $_[0]; # first parameter
    # can we open the file?
    unless (open(PROSITEFILE, $filename))
    {
        die "ERROR: Cannot open file ", $filename, "!\n";
    }

    while (<PROSITEFILE>)
    {
        chomp;
        $pattern .= $_;
    } # while

    close PROSITEFILE;

    # check if the pattern is not empty
    if ($pattern eq "")
    {
        die "ERROR: PROSITE file ", $filename, " contains empty pattern!\n";
    }
    return 0;
}

#===============================================================================
sub checkDataConsistency
#===============================================================================
{
    # check if the FASTA sequences relate to the NSEST data
    foreach my $nsest (@nsest_array)
    {
        my $name   = @$nsest[0];
        my $start  = @$nsest[1];
        my $end    = @$nsest[2];
        my $strand = @$nsest[3];
        my $type   = @$nsest[4];
        # run a three-step check if data is consistent
        # step 1. does the FASTA sequence exist for the NSEST line?
        if ($fasta_hash{$name} eq "")
        {
            die "ERROR : there is no FASTA sequence for NSEST\n",
                "Name  : ", $name,   "\n",
                "Start : ", $start,  "\n",
                "End   : ", $end,    "\n",
                "Strand: ", $strand, "\n",
                "Type  : ", $type,   "\n";
        }
        # step 2. is start > end?
        if ($start > $end)
        {
            die "ERROR : NSEST record contains corrupt data (start > end)\n",
                "Name  : ", $name,   "\n",
                "Start : ", $start,  "\n",
                "End   : ", $end,    "\n",
                "Strand: ", $strand, "\n",
                "Type  : ", $type,   "\n";
        }
        # step 3. is the length of an exon a multiple of 3?
        if (($type eq "exon") and ((($end - $start + 1) % 3) != 0))
        {
            die "ERROR: NSEST sequence ", $name, " starts at " , $start, " and finishes at ", $end, ", the length is ", ($end - $start + 1), " not multiple of 3\n";
        }
        # step 4. is the FASTA sequence long enough?
        if ($end > length($fasta_hash{$name}))
        {
            die "ERROR: FASTA sequence ", $name, " is too short (", length($fasta_hash{$name}), ") for the NSEST data (", $end, ")\n";
        }
    } # foreach
    return 0; # success
}

#===============================================================================
sub findExonsSequences
#===============================================================================
{
    # sort the NSEST array first
    @nsest_array = sort byNameType (@nsest_array);
    foreach my $nsest (@nsest_array)
    {
        my $type = @$nsest[4];
        if ($type eq "exon")
        {
            my $name          = @$nsest[0];
            my $start         = @$nsest[1];
            my $end           = @$nsest[2];
            my $strand        = @$nsest[3];
            my $orig_sequence = $fasta_hash{$name};
            #########################################################################
            # is this a reverse strand? There are combinations of three things to do:
            # reverse + complement + substring                        => !'s
            # reverse              + substring                        => !'s
            #           complement + substring                        => !'s
            #                        substring                        => !'s
            #                        substring + reverse              => OK...
            #                        substring           + complement => !'s
            #                        substring + reverse + complement => !'s
            # Should we try different frames then???
            if ($strand eq "R")
            {
                # make a reverse of the DNA sequence
                #$orig_sequence = reverse($orig_sequence);
                # substitute all bases by their complements (A->T, T->A, G->C, C->G)
                #$orig_sequence =~ tr/ACGT/TGCA/;
            }
            my $sub_sequence  = substr($orig_sequence, $start - 1, $end - $start + 1);
            if ($strand eq "R")
            {
                # make a reverse of the DNA sequence
                $sub_sequence = reverse($sub_sequence);
                # substitute all bases by their complements (A->T, T->A, G->C, C->G)
                #$sub_sequence =~ tr/ACGT/TGCA/;
            }
            #########################################################################
            @$nsest[5] = translateDNA2Protein($sub_sequence);
            print $name, "\t", $start, "\t", $end, "\t", $strand, "\t", $type, "\t", @$nsest[5], "\n";

        }
    } # foreach
    return 0; # success
}

#===============================================================================
sub findExonsMatchingPattern
#===============================================================================
{
    my $prosite_pattern = $_[0]; # first parameter
    my $perl_pattern = convertPROSITE2Perl($prosite_pattern);

    foreach my $nsest (@nsest_array)
    {
        my $type = @$nsest[4];
        if ($type eq "exon")
        {
            my $name             = @$nsest[0];
            my $start            = @$nsest[1];
            my $end              = @$nsest[2];
            my $strand           = @$nsest[3];
            my $search_sequence  = @$nsest[5]; # start the search from the beginning
            my $pattern_start    = 0;          # assuming there is no pattern
            my $pattern_end      = 0;          # assuming there is no pattern
            my $search_offset    = 0;          # initial value of relative coordinate
            my $no_more_found    = 0;          # initial value is false
            do
            {
                if ($search_sequence =~ /$perl_pattern/)
                {
                    # we found something but what exactly...
                    my $matching = $&; # sequence that matches the pattern
                    $pattern_start = index($search_sequence, $matching); # zero-based
                    if ($pattern_start >= 0)
                    {
                        $pattern_end = $pattern_start + length ($matching) - 1; # zero-based

                        # fill the result array
                        my @result_line;
                        $result_line[0] = @$nsest[0];                      # name
                        $result_line[1] = $search_offset + $pattern_start; # start of the matching sequence
                        $result_line[2] = $search_offset + $pattern_end;   # end of the matching sequence
                        $result_line[3] = $matching;                       # proteine that matches the pattern
                        push @result_array, \@result_line;

                        # change the starting point for the search
                        # start the search again after the first character of matched string
                        $search_offset += $pattern_start + 1;
                        # reduce the sequence for the next search
                        $search_sequence = substr ($search_sequence, $pattern_start + 1);
                    } # if pattern actually located
                } # if match found
                else
                {
                    # there are no more occurences of pattern
                    $no_more_found = 1;
                } # else
            } # do
            until ($no_more_found)
        } # if this is exon
    } # foreach
    return 0; # success
}

#===============================================================================
sub byNameType
#===============================================================================
{
    return (($a->[0] cmp $b->[0]) or  # name
            ($a->[4] cmp $b->[4]) or  # type
            ($a->[1] <=> $b->[1]) or  # start
            ($a->[2] <=> $b->[2]));   # end
}


#===============================================================================
sub translateDNA2Protein
#===============================================================================
{
    my $dna_seq = $_[0]; # first parameter

    my %codes = ('AAA'=>'K', 'AAC'=>'N', 'AAG'=>'K', 'AAT'=>'N',
                 'ACA'=>'T', 'ACC'=>'T', 'ACG'=>'T', 'ACT'=>'T',
                 'AGA'=>'R', 'AGC'=>'S', 'AGG'=>'R', 'AGT'=>'S',
                 'ATA'=>'I', 'ATC'=>'I', 'ATG'=>'M', 'ATT'=>'I',
                 'CAA'=>'Q', 'CAC'=>'H', 'CAG'=>'Q', 'CAT'=>'H',
                 'CCA'=>'P', 'CCC'=>'P', 'CCG'=>'P', 'CCT'=>'P',
                 'CGA'=>'R', 'CGC'=>'R', 'CGG'=>'R', 'CGT'=>'R',
                 'CTA'=>'L', 'CTC'=>'L', 'CTG'=>'L', 'CTT'=>'L',
                 'GAA'=>'E', 'GAC'=>'D', 'GAG'=>'E', 'GAT'=>'D',
                 'GCA'=>'A', 'GCC'=>'A', 'GCG'=>'A', 'GCT'=>'A',
                 'GGA'=>'G', 'GGC'=>'G', 'GGG'=>'G', 'GGT'=>'G',
                 'GTA'=>'V', 'GTC'=>'V', 'GTG'=>'V', 'GTT'=>'V',
                 'TAA'=>'!', 'TAC'=>'Y', 'TAG'=>'!', 'TAT'=>'Y',
                 'TCA'=>'S', 'TCC'=>'S', 'TCG'=>'S', 'TCT'=>'S',
                 'TGA'=>'!', 'TGC'=>'C', 'TGG'=>'W', 'TGT'=>'C',
                 'TTA'=>'L', 'TTC'=>'F', 'TTG'=>'L', 'TTT'=>'F');

    my $dna_len = length($dna_seq);

    if (($dna_len) % 3 != 0) {print "ERROR: Length of DNA sequence '", $dna_seq, "' is not multiple of 3\n"; exit}

    my $protein_seq = "";

    for (my $i = 0; $i < $dna_len; $i += 3)
    {
        my $codon = substr($dna_seq, $i, 3);
        if (! defined ($codes{$codon})) {print "ERROR: Illigal codon ", $codon, " in sequence ", $dna_seq, "\n"; exit;}
        $protein_seq .= $codes{$codon};
    } # for
    return $protein_seq;
}

#===============================================================================
sub convertPROSITE2Perl
#===============================================================================
{
    my $prosite_pattern = $_[0]; # first parameter

    my $perl_pattern = $prosite_pattern;
    # remove '-'
    $perl_pattern =~ s/-//g;
    # replace '(' with '{'
    $perl_pattern =~ s/\(/{/g;
    # replace ')' with '}'
    $perl_pattern =~ s/\)/}/g;
    return $perl_pattern;
}

#===============================================================================
sub htmlPrintHeader
#===============================================================================
{
    my $title = $_[0]; # HTML page title is the first parameter
    print "<HTML>\n";
    print "  <HEAD>\n";
    print "    <TITLE>", $title, "</TITLE>\n";
    print "  </HEAD>\n";
    print "  <BODY>\n";
    return 0;
}

#===============================================================================
sub htmlPrintFooter
#===============================================================================
{
    print "  </BODY>\n";
    print "</HTML>\n";
    return 0;
}

#===============================================================================
sub htmlPrintParagraph
#===============================================================================
{
    my $text = $_[0]; # text is the first parameter
    print "  <P>", $text, "\n";
    return 0;
}


#===============================================================================
sub htmlPrintTable
#===============================================================================
{
    my @table = @_;
    my $display_header = 1; # first row contains headers for the columns
    print "<TABLE BORDER='1'>\n";
    print "<CAPTION ALIGN='TOP'>Results</CAPTION>\n";
    foreach my $row (@table)
    {
        print "<TR>\n";
        foreach my $cell (@$row)
        {
            if ($display_header == 1) {print "<TH>";} else {print "<TD>";}
            print $cell;
        }
        $display_header = 0; # display column headers once only
        print "\n";
    }
    print "</TABLE>\n";
    return 0;
}


