#!/usr/bin/perl

use strict;
use warnings;

# Input and output file names
my $input_file = "specs_ms.mgf";
my $output_file = "specs_ms_sel41.mgf";
my $scan_numbers_file = "scan_numbers_sel41.txt";

# Read scan numbers from the file
open my $scan_numbers_fh, '<', $scan_numbers_file or die "Cannot open scan numbers file: $!";
my @scan_numbers = map { chomp; $_ } <$scan_numbers_fh>;
close $scan_numbers_fh;

# Open input and output files
open my $input_fh, '<', $input_file or die "Cannot open input file: $!";
open my $output_fh, '>', $output_file or die "Cannot open output file: $!";

# Variables to track whether to print the current spectrum
my $in_spectrum = 0;
my $current_scan_number;

# Process each line in the input file
while (<$input_fh>) {
    if (/^BEGIN IONS/) {
        $in_spectrum = 0;
    }
    elsif (/^SCANS=(\d+)/) {
        $current_scan_number = $1;
        $in_spectrum = 1 if grep { $_ eq $current_scan_number } @scan_numbers;
    }

    print $output_fh $_ if $in_spectrum;

    $in_spectrum = 0 if /^END IONS/;
}

# Close the files
close $input_fh;
close $output_fh;

print "Extraction completed. Output written to $output_file\n";
