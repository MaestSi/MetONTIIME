#!/usr/bin/perl

#
# Copyright 2019 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

use strict;
use warnings;

#
# READ LENGTH THRESHOLD
# filter out reads with length > READ_LENGTH_MAX or length < READ_LENGTH_MIN
#

my ($READ_LENGTH_MIN, $READ_LENGTH_MAX) = @ARGV;


while (<STDIN>)
{
    chomp;
    my $name = $_;
    my $seq = <STDIN>;
    chomp $seq;
    my $seq_len = length($seq);
    my $symb = <STDIN>;
    chomp $symb;
    my $qual = <STDIN>;
    chomp $qual;
    

    #
    # READ LENGTH THRESHOLD
    #
    next if ( $seq_len > $READ_LENGTH_MAX || $seq_len < $READ_LENGTH_MIN );

    print join "\n", $name, $seq, $symb, $qual, '';
}
