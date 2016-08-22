#!/usr/bin/perl -w

###############################################################################
#
# This Perl script implements the various tokenization heuristics studied in
# [1]. These heuristics are specifically designed for ad-hoc information
# retrieval from biomedical text, and we found that for different types of 
# queries, different heuristics should be used to achieve the best performance.
# The details are discussed in [1].
# 
# The script can be executed from the command line as follows:
#     $ BioTokenizer.pl -i input -o output [-h] [-t <S|V>] [-b <0|1|2|3>] [-n <h|s|j>] [-g] [-s <p|l|s>]
# 
# When the option "-h" is chosen, some help information will be printed to the 
# screen.
# 
# Use the options "-i" and "-o" to specify the input and the output files, 
# respectively. There two files must be specified. In mose cases, you can just
# specify these two files, and ignore the other settings. The default 
# tokenization strategies will then be used.
# 
# Use the option "-t" to specify the query type. If your queries consiste of 
# only gene/protein symbols, set "-t" to "S" (for symbolic queries). If your 
# queries contain full protein/gene names and other English words, set "-t" to
# "V" (for verbose queries). Note that setting "-t" to "V" is the same as 
# using the default settings.
#
# You can try more flexible combinations of tokenization strategies by using
# the "-b", "-n", "-g" and "-s" options. Use the option "-b" to specify the
# set of break points (see [1] for details). When "-b" is set to "0", no break
# point is considered; otherwise, "1", "2" and "3" correspond to BP1, BP2 and 
# BP3 in [1]. Use the option "-n" to specify the break point normalization 
# method. Use the option "-g" to specify whether or not Greek alphabet
# normalization should be performed. Use the option "-s" to specify the 
# stemming method. "p" stands for the Porter stemmer, "l" stands for the Lovins
# stemmer, and "s" stands for the S stemmer.
#
# The input file is a plain text file. To make the script compatible with 
# Lemur/Indri document files, the following lines are allowed in the input 
# file, and will not be tokenized.
#
# <DOC>
# <DOCNO>xxxx</DOCNO>
# </DOC>
#
# "xxxx" is any string, and is used as a document ID.
# 
# With the current implementation of this script, it takes about 1 minute to 
# process 5M of text on a machine with a 3.0GHz CPU and 1G of memory.
#
# The code for the Porter stemming algorithm and the Lovins stemming algorithm
# is taken from the implementation by the following authors:
#
# Porter: 
#     Martin Porter
#     http://www.tartarus.org/~martin/PorterStemmer/perl.txt
# 
# Lovins:
#     Lovins-stemmer.pl Version 1.1
#     Copyright 2001 Gordon W. Paynter (paynter@cs.waikato.ac.nz).
#     Distributed under the GNU General Public License (Version 2).
###############################################################################

use strict;

use Getopt::Std;

my(%opt) = ();
getopts('ht:i:o:b:n:gs:', \%opt);

my($usage) = "Usage: $0 -i input -o output [-h] [-t <S|V>] [-b <0|1|2|3>] [-n <h|s|j>] [-g] [-s <p|l|s>]\n"
    ."\t-i: input\n"
    ."\t-o: output\n"
    ."\t-h: help information\n"
    ."\t-t: query type [S|V]\n"
    ."\t-b: break point set [0|1|2|3]\n"
    ."\t-n: break point normalization method [h|s|j]\n"
    ."\t-g: Greek alphabet normalization [on|off]\n"
    ."\t-s: stemming method [p(orter)|l(ovins)|s]\n";

if (defined($opt{h})) {
    print STDERR $usage;
    exit;
}

if (!defined($opt{i}) || !defined($opt{o})) {
    print STDERR $usage;
    exit;
}

my($bp, $norm, $grk, $stem);

if (defined($opt{t})) {

    if ($opt{t} eq "S") {
	$bp = 1;
	$norm = "j";
	$grk = 1;
	$stem = "";
    } elsif ($opt{t} eq "V") {
	$bp = 1;
	$norm = "s";
	$grk = 0;
	$stem = "p";
    } else {
	print STDERR "The query type must be S(Symbolic) or V(Verbose)!\n";
	print STDERR $usage;
	exit;
    }

} else {

    if (defined($opt{b})) {
	if ($opt{b} ne "0" && $opt{b} ne "1" && $opt{b} ne "2" 
	    && $opt{b} ne "3") {
	    print STDERR "The break point set must be 0, 1, 2 or 3!\n";
	    print STDERR $usage;
	    exit;
	} else {
	    $bp = $opt{b};
	}

	if ($opt{b} eq "0") {
	    $norm = "";
	    $grk = 0;
	    $stem = "";
	} else {
	    if (!defined($opt{n})) {
		print STDERR "If the break point set is not 0, a normalization method must be specified!\n";
		print STDERR $usage;
		exit;
	    }
	    
	    if ($opt{n} ne "h" && $opt{n} ne "s" && $opt{n} ne "j") {
		print STDERR "The normalization method must be 'h', 's' or 'j'!\n";
		print STDERR $usage;
		exit;
	    }

	    $norm = $opt{n};

	    if (defined($opt{g})) {
		$grk = 1;
	    } else {
		$grk = 0;
	    }

	    if (defined($opt{s})) {
		if ($opt{s} eq "p" || $opt{s} eq "l"
		    || $opt{s} eq "s") {
		    $stem = $opt{s};
		} else {
		    print STDERR "The stemming method must be 'p', 'l' or 's'!\n";
		    print STDERR $usage;
		    exit;
		}
	    } else {
		$stem = "";
	    }
	}	
	
    } else {
	$bp = 1;
	$norm = "s";
	$grk = 0;
	$stem = "p";	
    }
}

my($i);

my(%greek) = ("alpha" => "a",
	      "beta" => "b",
	      "gamma" => "g",
	      "delta" => "d",
	      "epsilon" => "e",
	      "zeta" => "z",
	      "eta" => "e",
	      "theta" => "th",
	      "iota" => "i",
	      "kappa" => "k",
	      "lambda" => "l",
	      "mu" => "m",
	      "nu" => "n",
	      "xi" => "x",
	      "omicron" => "o",
	      "pi" => "p",
	      "rho" => "r",
	      "sigma" => "s",
	      "tau" => "t",
	      "upsilon" => "u",
	      "phi" => "ph",
	      "chi" => "ch",
	      "psi" => "ps",
	      "omega" => "o");
	      
my(%step2list) = ();
my(%step3list) = ();
my($c, $v, $C, $V, $mgr0, $meq1, $mgr1, $_v);
my(%ending_list) = ();

if ($stem ne "") {
    initialise();
}

open(IN, "</home/akshay/2630847");
open(OUT, ">/home/akshay/output");
my($line);
while($line = <IN>) {
    if ($line =~ /^<DOCNO/ || $line =~ /^<DOC>/
	|| $line =~ /^<\/DOC>/) {
	print OUT $line;
	next;
    }
    if ($line =~ /^<TITLE/ ||$line =~ /^<\/TITLE>/) {
	print OUT $line;
	next;
    }
    if ($line =~ /^<TEXT/ ||$line =~ /^<\/TEXT>/) {
	print OUT $line;
	next;
    }
    chomp($line);
    
    $line = " ".$line." ";
    
    # ============================================
    # removal of non-functional characters
    
    # heuristic rule 1
    $line =~ s/[\!\"\#\$\%\&\*\<\=\>\?\@\\\|\~]/ /g;   
    # heuristic rule 2
    $line =~ s/[\.\:\;\,] / /g;
    # heuristic rule 3
    $line =~ s/ \(([^\)]*)\) / $1 /g;
    $line =~ s/ \[([^\)]*)\] / $1 /g;
    # repeat heuristic rule 3 for nested parentheses
    $line =~ s/ \(([^\)]*)\) / $1 /g;
    $line =~ s/ \[([^\)]*)\] / $1 /g;
    # heuristic rule 4
    $line =~ s/ [\']/ /g;
    $line =~ s/[\`] / /g;
    # heuristic rule 5
    $line =~ s/[\'][st] / /g;
    # repeat heuristic rule 2 for cases where these symbols were originally followed by other symbols such as parentheses
    $line =~ s/[\.\:\;\,] / /g;
    # heuristic rule 6
    $line =~ s/[\/]+ / /g;
    # ============================================

    $line =~ s/^\s+//;
    $line =~ s/\s+$//;
    
    my(@tokens) = split(/\s+/, $line);
    my($token);

    foreach $token (@tokens) {
	# ============================================
	# split the token by break points
	
	my(@subtokens) = ();

	if ($bp eq "0") {
	    push(@subtokens, $token);
	} elsif ($bp eq "1") {
	    while($token =~ /[^\(\)\[\]\-\_\/]+/) {
		push(@subtokens, $&);
		$token = "$'";
	    }
	} elsif ($bp eq "2") {
	    while($token =~ /[A-Za-z0-9]+/) {
		push(@subtokens, $&);
		$token = "$'";
	    }
	} else {
	    while($token =~ /([A-Z][a-z]+)|([A-Z]+)|([a-z]+)|([0-9]+)/) {
		push(@subtokens, $&);
		$token = "$'";
	    }
	}

	# ============================================
	# transform to lowercase
	for($i = 0; $i < scalar(@subtokens); $i++) {
	    $subtokens[$i] = lc($subtokens[$i]);
	}	

	# ============================================
	# we do the Greek alphabet normalization before normalizing the break points
	if ($grk) {
	    for($i = 0; $i < scalar(@subtokens); $i++) {
		my($rest) = $subtokens[$i];
		my($pre) = "";
		while($rest =~ /[a-z]+/) {
		    $rest = "$'";
		    $pre .= "$`";
		    my($substring) = $&;
		    if (defined($greek{$substring})) {
			$pre .= $greek{$substring};
		    } else {
			$pre .= $substring;
		    }
		}
		$subtokens[$i] = $pre.$rest;
	    }
	}

	# ============================================
	# normalize the break points
	
	if ($bp eq "0") {
	    $token = $subtokens[0];
	} elsif ($norm eq "h") {
	    $token = join('-', @subtokens);
	} elsif ($norm eq "s") {

	    # ============================================
	    # stemming
	    if ($stem ne "") {
		for($i = 0; $i < scalar(@subtokens); $i++) {
		    if ($stem eq "p") {
			$subtokens[$i] = PorterStem($subtokens[$i]);
		    } elsif ($stem eq "l") {
			$subtokens[$i] = LovinsStem($subtokens[$i]);
		    } else {
			$subtokens[$i] = SStem($subtokens[$i]);
		    }
		}
	    }

	    $token = join(' ', @subtokens);
	} else {
	    $token = join('', @subtokens);
	}

	# ============================================
	# stemming
	
	if ($stem ne "") {
	    if ($bp eq "0" || 
		($bp ne "0" && ($norm eq "h" || $norm eq "j"))) {
		if ($stem eq "p") {
		    $token = PorterStem($token);
		} elsif ($stem eq "l") {
		    $token = LovinsStem($token);
		} else {
		    $token = SStem($token);
		}
	    }
	}
	
	# ============================================

	print OUT $token, " ";
    }

    print OUT "\n";
}

close(IN);
close(OUT);

sub PorterStem {
    my ($stem, $suffix, $firstch);
    my $w = shift;
    if (length($w) < 3) { return $w; } # length at least 3
    # now map initial y to Y so that the patterns never treat it as vowel:
    $w =~ /^./; $firstch = $&;
    if ($firstch =~ /^y/) { $w = ucfirst $w; }
   
    # Step 1a
    if ($w =~ /(ss|i)es$/) { $w=$`.$1; }
    elsif ($w =~ /([^s])s$/) { $w=$`.$1; }
    # Step 1b
    if ($w =~ /eed$/) { if ($` =~ /$mgr0/o) { chop($w); } }
    elsif ($w =~ /(ed|ing)$/) {
	$stem = $`;
	if ($stem =~ /$_v/o) {
	    $w = $stem;
	    if ($w =~ /(at|bl|iz)$/) { $w .= "e"; }
	    elsif ($w =~ /([^aeiouylsz])\1$/) { chop($w); }
	    elsif ($w =~ /^${C}${v}[^aeiouwxy]$/o) { $w .= "e"; }
	}
    }
    # Step 1c
    if ($w =~ /y$/) { $stem = $`; if ($stem =~ /$_v/o) { $w = $stem."i"; } }
    
    # Step 2
    if ($w =~ /(ational|tional|enci|anci|izer|bli|alli|entli|eli|ousli|ization|ation|ator|alism|iveness|fulness|ousness|aliti|iviti|biliti|logi)$/) {
	$stem = $`; $suffix = $1;
	if ($stem =~ /$mgr0/o) { $w = $stem . $step2list{$suffix}; }
    }
   
    # Step 3
    
    if ($w =~ /(icate|ative|alize|iciti|ical|ful|ness)$/) {
	$stem = $`; $suffix = $1;
	if ($stem =~ /$mgr0/o) { $w = $stem . $step3list{$suffix}; }
    }
   
    # Step 4
    
    if ($w =~ /(al|ance|ence|er|ic|able|ible|ant|ement|ment|ent|ou|ism|ate|iti|ous|ive|ize)$/) {
	$stem = $`; if ($stem =~ /$mgr1/o) { $w = $stem; }
    } elsif ($w =~ /(s|t)(ion)$/) {
	$stem = $` . $1; if ($stem =~ /$mgr1/o) { $w = $stem; }
    }
        
    #  Step 5
   
    if ($w =~ /e$/) {
	$stem = $`;
	if ($stem =~ /$mgr1/o or
	    ($stem =~ /$meq1/o and not $stem =~ /^${C}${v}[^aeiouwxy]$/o))
	{ $w = $stem; }
    }
    if ($w =~ /ll$/ and $w =~ /$mgr1/o) { chop($w); }
   
    # and turn initial Y back to y
    if ($firstch =~ /^y/) { $w = lcfirst $w; }
    return $w;
}

sub LovinsStem {
    my $word = shift @_;

    # the word length
    my $wlen = length $word;
    return $word unless ($wlen > 2);

    # Phase 1: search for a match between the phrase and the list of endings
    my ($prefix, $suffix, $prelen, $suflen);

    # The longest suffix is 11 characters long
    if ($wlen <= 13) {
	$prelen = 2;
	$prefix = substr $word, 0, 2;
	$suffix = substr $word, 2;
	$suflen = length $suffix;
    } else {
	$prelen = $wlen - 11;
	$prefix = substr $word, 0, $prelen;
	$suffix = substr $word, $prelen;
	$suflen = length $suffix;
    }

    my ($condition_code, $stem);
    for (;;) {
	if ($condition_code = $ending_list{$suffix}) {
	    # Stem if context-sensitive rules are satisfied
	    if (($condition_code eq 'A')  
		# I've sorted these codes in the order of the frequency of their
	        # in the brown corpus (more or less).
		|| (($condition_code eq 'B') && ($prelen >= 3))
		|| (($condition_code eq 'W') && ($prefix !~ /[su]$/o))
		|| (($condition_code eq 'E') && ($prefix !~ /e$/o))

		|| (($condition_code eq 'N') && ($prefix =~ /([^s]..|.s..)$/o))
		|| (($condition_code eq 'F') && ($prelen >= 3) && ($prefix !~ /e$/o))

		|| (($condition_code eq 'Q') && ($prelen >= 3) && ($prefix !~ /[ln]$/o))
		|| (($condition_code eq 'C') && ($prelen >= 4))
		|| (($condition_code eq 'BB') && ($prelen >= 3) && ($prefix !~ /(met|ryst)$/o))

		|| (($condition_code eq 'S') && ($prefix =~ /(dr|[^t]t)$/o))
		|| (($condition_code eq 'T') && ($prefix =~ /(s|[^o]t)$/o))
		|| (($condition_code eq 'X') && ($prefix =~ /(l|i|u.e)$/o))
		|| (($condition_code eq 'I') && ($prefix !~ /[oe]$/o))
		|| (($condition_code eq 'P') && ($prefix !~ /c$/o))

		|| (($condition_code eq 'M') && ($prefix !~ /[aecm]$/o))
		|| (($condition_code eq 'L') && ($prefix !~ /(u|x|[^o]s)$/o))
		|| (($condition_code eq 'O') && ($prefix =~ /[li]$/o))
		|| (($condition_code eq 'AA') && ($prefix =~ /(d|f|ph|th|l|er|or|es|t)$/o))
		|| (($condition_code eq 'V') && ($prefix =~ /c$/o))

		|| (($condition_code eq 'R') && ($prefix =~ /[nr]$/o))
		|| (($condition_code eq 'Y') && ($prefix =~ /in$/o))
		|| (($condition_code eq 'G') && ($prelen >= 3) && ($prefix =~ /f$/o))
		|| (($condition_code eq 'K') && ($prelen >= 3) && ($prefix =~ /(l|i|u.e)$/o))
		|| (($condition_code eq 'U') && ($prefix =~ /[lmnr]$/o))
		|| (($condition_code eq 'D') && ($prelen >= 5))
		|| (($condition_code eq 'H') && ($prefix =~ /(t|ll)$/o))
		|| (($condition_code eq 'J') && ($prefix !~ /[ae]$/o))
		|| (($condition_code eq 'Z') && ($prefix !~ /f$/o))
		|| (($condition_code eq 'CC') && ($prefix =~ /l$/o)) ) {
		
		$word = $prefix;
		$wlen = $prelen;
		last;
	    }
	}
	
	# Try the next shorter suffix
	$prefix .= substr $suffix, 0, 1;
	$prelen++;
	$suffix = substr $suffix, 1;
	$suflen--;
	last if ($suflen == 0);

       
    }

    # Phase 2: recoding suffixes
    my $last1 = substr $word, ($wlen - 1);
    if ($last1 eq 't') {
	$word =~ s/tt$/t/o;
	$word =~ s/uct$/uc/o;
	$word =~ s/umpt$/um/o;
	$word =~ s/rpt$/rb/o;
	$word =~ s/mit$/mis/o;
	$word =~ s/ert$/ers/o;
	$word =~ s/^et$/es/o;
	$word =~ s/([^n])et$/$1es/o;
        $word =~ s/yt$/ys/o;
    } elsif ($last1 eq 'r') {
	$word =~ s/rr$/r/o;
	$word =~ s/istr$/ister/o;
	$word =~ s/metr$/meter/o;
	$word =~ s/^her$/hes/o;
	$word =~ s/([^pt])her$/$1hes/o;
    } elsif ($last1 eq 'd') {
	$word =~ s/dd$/d/o;
	$word =~ s/uad$/uas/o;
	$word =~ s/vad$/vas/o;
	$word =~ s/cid$/cis/o;
	$word =~ s/lid$/lis/o;
	$word =~ s/erid$/eris/o;
	$word =~ s/pand$/pans/o;
	$word =~ s/^end$/ens/o;
	$word =~ s/([^s])end$/$1ens/o;
	$word =~ s/ond$/ons/o;
	$word =~ s/lud$/lus/o;
	$word =~ s/rud$/rus/o;
	$word =~ s/^end$/ens/o;
	$word =~ s/([^m])end$/$1ens/o;
    } elsif ($last1 eq 'n') {
	$word =~ s/nn$/n/o;
    } elsif ($last1 eq 'l') {
	$word =~ s/ll$/l/o;
	
	# J. Jiang
	# $word =~ s/^ul$/$1l/o;

	$word =~ s/([^aio])ul$/$1l/o;
    } elsif ($last1 eq 'm') {
	$word =~ s/mm$/m/o;
    } elsif ($last1 eq 's') {
	$word =~ s/ss$/s/o;
	$word =~ s/urs$/ur/o;
    } elsif ($last1 eq 'g') {
	$word =~ s/gg$/g/o;
    } elsif ($last1 eq 'v') {
	$word =~ s/iev$/ief/o;
	$word =~ s/olv$/olut/o;
    } elsif ($last1 eq 'p') {
	$word =~ s/pp$/p/o;
    } elsif ($last1 eq 'b') {
	$word =~ s/bb$/b/o;
    } elsif ($last1 eq 'x') {
	$word =~ s/bex$/bic/o;
	$word =~ s/dex$/dic/o;
	$word =~ s/pex$/pic/o;
	$word =~ s/tex$/tic/o;
	$word =~ s/ax$/ac/o;
	$word =~ s/ex$/ec/o;
	$word =~ s/ix$/ic/o;
	$word =~ s/lux$/luc/o;
    } elsif ($last1 eq 'z') {
	$word =~ s/yz$/ys/o;
    }

    # output
    return $word;
}

sub SStem{
    my($word) = $_[0];
    if ($word =~ /ies$/) {
	$word =~ s/^ies$/y/;
	$word =~ s/([^ea])ies/$1y/;	
    }
    elsif ($word =~ /es$/) {
	$word =~ s/^es/e/;
	$word =~ s/([^aeo])es/$1e/;	
    }
    elsif ($word =~ /s$/) {
	$word =~ s/([^us])s$/$1/;
    }

    return $word;
}

sub initialise {

    %step2list =
	( 'ational'=>'ate', 'tional'=>'tion', 'enci'=>'ence', 'anci'=>'ance', 'izer'=>'ize', 'bli'=>'ble',
	  'alli'=>'al', 'entli'=>'ent', 'eli'=>'e', 'ousli'=>'ous', 'ization'=>'ize', 'ation'=>'ate',
	  'ator'=>'ate', 'alism'=>'al', 'iveness'=>'ive', 'fulness'=>'ful', 'ousness'=>'ous', 'aliti'=>'al',
	  'iviti'=>'ive', 'biliti'=>'ble', 'logi'=>'log');
    
    %step3list =
	('icate'=>'ic', 'ative'=>'', 'alize'=>'al', 'iciti'=>'ic', 'ical'=>'ic', 'ful'=>'', 'ness'=>'');
    
    
    $c =    "[^aeiou]";          # consonant
    $v =    "[aeiouy]";          # vowel

    # #####################################################
    # modified by J. Jiang

    # $C =    "${c}[^aeiouy]*";    # consonant sequence
    $C =    "$c"."[^aeiouy]*";    # consonant sequence
    
    # $V =    "${v}[aeiou]*";      # vowel sequence
    $V =    "$v"."[aeiou]*";      # vowel sequence
    # #####################################################
    
    $mgr0 = "^(${C})?${V}${C}";               # [C]VC... is m>0
    $meq1 = "^(${C})?${V}${C}(${V})?" . '$';  # [C]VC[V] is m=1
    $mgr1 = "^(${C})?${V}${C}${V}${C}";       # [C]VCVC... is m>1
    $_v   = "^(${C})?${v}";                   # vowel in stem

    %ending_list = (
		    'a' => 'A',
		    
		    'ae' => 'A',
		    'al' => 'BB',
		    'ar' => 'X',
		    'as' => 'B',

		    'acy' => 'A',
		    'age' => 'B',
		    'aic' => 'A',
		    'als' => 'BB',
		    'ant' => 'B',
		    'ars' => 'O',
		    'ary' => 'F',
		    'ata' => 'A',
		    'ate' => 'A',

		    'able' => 'A',
		    'ably' => 'A',
		    'ages' => 'B',
		    'ally' => 'B',
		    'ance' => 'B',
		    'ancy' => 'B',
		    'ants' => 'B',
		    'aric' => 'A',
		    'arly' => 'K',
		    'ated' => 'I',
		    'ates' => 'A',
		    'atic' => 'B',
		    'ator' => 'A',

		    'acies' => 'A',
		    'acity' => 'A',
		    'aging' => 'B',
		    'aical' => 'A',
		    'alist' => 'A',
		    'alism' => 'B',
		    'ality' => 'A',
		    'alize' => 'A',
		    'allic' => 'BB',
		    'anced' => 'B',
		    'ances' => 'B',
		    'antic' => 'C',
		    'arial' => 'A',
		    'aries' => 'A',
		    'arily' => 'A',
		    'arity' => 'B',
		    'arize' => 'A',
		    'aroid' => 'A',
		    'ately' => 'A',
		    'ating' => 'I',
		    'ation' => 'B',
		    'ative' => 'A',
		    'ators' => 'A',
		    'atory' => 'A',
		    'ature' => 'E',

		    'aceous' => 'A',
		    'acious' => 'B',
		    'action' => 'G',
		    'alness' => 'A',
		    'ancial' => 'A',
		    'ancies' => 'A',
		    'ancing' => 'B',
		    'ariser' => 'A',
		    'arized' => 'A',
		    'arizer' => 'A',
		    'atable' => 'A',
		    'ations' => 'B',
		    'atives' => 'A',
		    
		    'ability' => 'A',
		    'aically' => 'A',
		    'alistic' => 'B',
		    'alities' => 'A',
		    'ariness' => 'E',
		    'aristic' => 'A',
		    'arizing' => 'A',
		    'ateness' => 'A',
		    'atingly' => 'A',
		    'ational' => 'B',
		    'atively' => 'A',
		    'ativism' => 'A',
		    
		    'ableness' => 'A',
		    'arizable' => 'A',
		    
		    'allically' => 'C',
		    'antaneous' => 'A',
		    'antiality' => 'A',
		    'arisation' => 'A',
		    'arization' => 'A',
		    'ationally' => 'B',
		    'ativeness' => 'A',
		    
		    'antialness' => 'A',
		    'arisations' => 'A',
		    'arizations' => 'A',

		    'alistically' => 'B',
		    'arizability' => 'A',

		    'e' => 'A',
		    
		    'ed' => 'E',
		    'en' => 'F',
		    'es' => 'E',
		    
		    'eal' => 'Y',
		    'ear' => 'Y',
		    'ely' => 'E',
		    'ene' => 'E',
		    'ent' => 'C',
		    'ery' => 'E',
		    'ese' => 'A',
		    
		    'ealy' => 'Y',
		    'edly' => 'E',
		    'eful' => 'A',
		    'eity' => 'A',
		    'ence' => 'A',
		    'ency' => 'A',
		    'ened' => 'E',
		    'enly' => 'E',
		    'eous' => 'A',
		    
		    'early' => 'Y',
		    'ehood' => 'A',
		    'eless' => 'A',
		    'elily' => 'A',
		    'ement' => 'A',
		    'enced' => 'A',
		    'ences' => 'A',
		    'eness' => 'E',
		    'ening' => 'E',
		    'ental' => 'A',
		    'ented' => 'C',
		    'ently' => 'A',
		    
		    'eature' => 'Z',
		    'efully' => 'A',
		    'encies' => 'A',
		    'encing' => 'A',
		    'ential' => 'A',
		    'enting' => 'C',
		    'entist' => 'A',
		    'eously' => 'A',

		    'elihood' => 'E',
		    'encible' => 'A',
		    'entally' => 'A',
		    'entials' => 'A',
		    'entiate' => 'A',
		    'entness' => 'A',

		    'entation' => 'A',
		    'entially' => 'A',
		    'eousness' => 'A',

		    'eableness' => 'E',
		    'entations' => 'A',
		    'entiality' => 'A',
		    'entialize' => 'A',
		    'entiation' => 'A',

		    'entialness' => 'A',

		    'ful' => 'A',

		    'fully' => 'A',

		    'fulness' => 'A',

		    'hood' => 'A',

		    'i' => 'A',

		    'ia' => 'A',
		    'ic' => 'A',
		    'is' => 'A',

		    'ial' => 'A',
		    'ian' => 'A',
		    'ics' => 'A',
		    'ide' => 'L',
		    'ied' => 'A',
		    'ier' => 'A',
		    'ies' => 'P',
		    'ily' => 'A',
		    'ine' => 'M',
		    'ing' => 'N',
		    'ion' => 'Q',
		    'ish' => 'C',
		    'ism' => 'B',
		    'ist' => 'A',
		    'ite' => 'AA',
		    'ity' => 'A',
		    'ium' => 'A',
		    'ive' => 'A',
		    'ize' => 'F',

		    'ials' => 'A',
		    'ians' => 'A',
		    'ible' => 'A',
		    'ibly' => 'A',
		    'ical' => 'A',
		    'ides' => 'L',
		    'iers' => 'A',
		    'iful' => 'A',
		    'ines' => 'M',
		    'ings' => 'N',
		    'ions' => 'B',
		    'ious' => 'A',
		    'isms' => 'B',
		    'ists' => 'A',
		    'itic' => 'H',
		    'ized' => 'F',
		    'izer' => 'F',

		    'ially' => 'A',
		    'icant' => 'A',
		    'ician' => 'A',
		    'icide' => 'A',
		    'icism' => 'A',
		    'icist' => 'A',
		    'icity' => 'A',
		    'idine' => 'I',
		    'iedly' => 'A',
		    'ihood' => 'A',
		    'inate' => 'A',
		    'iness' => 'A',
		    'ingly' => 'B',
		    'inism' => 'J',
		    'inity' => 'CC',
		    'ional' => 'A',
		    'ioned' => 'A',
		    'ished' => 'A',
		    'istic' => 'A',
		    'ities' => 'A',
		    'itous' => 'A',
		    'ively' => 'A',
		    'ivity' => 'A',
		    'izers' => 'F',
		    'izing' => 'F',

		    'ialist' => 'A',
		    'iality' => 'A',
		    'ialize' => 'A',
		    'ically' => 'A',
		    'icance' => 'A',
		    'icians' => 'A',
		    'icists' => 'A',
		    'ifully' => 'A',
		    'ionals' => 'A',
		    'ionate' => 'D',
		    'ioning' => 'A',
		    'ionist' => 'A',
		    'iously' => 'A',
		    'istics' => 'A',
		    'izable' => 'E',

		    'ibility' => 'A',
		    'icalism' => 'A',
		    'icalist' => 'A',
		    'icality' => 'A',
		    'icalize' => 'A',
		    'ication' => 'G',
		    'icianry' => 'A',
		    'ination' => 'A',
		    'ingness' => 'A',
		    'ionally' => 'A',
		    'isation' => 'A',
		    'ishness' => 'A',
		    'istical' => 'A',
		    'iteness' => 'A',
		    'iveness' => 'A',
		    'ivistic' => 'A',
		    'ivities' => 'A',
		    'ization' => 'F',
		    'izement' => 'A',

		    'ibleness' => 'A',
		    'icalness' => 'A',
		    'ionalism' => 'A',
		    'ionality' => 'A',
		    'ionalize' => 'A',
		    'iousness' => 'A',
		    'izations' => 'A',

		    'ionalness' => 'A',
		    'istically' => 'A',
		    'itousness' => 'A',
		    'izability' => 'A',
		    'izational' => 'A',

		    'izationally' => 'B',

		    'ly' => 'B',

		    'less' => 'A',
		    'lily' => 'A',

		    'lessly' => 'A',

		    'lessness' => 'A',

		    'ness' => 'A',

		    'nesses' => 'A',

		    'o' => 'A',

		    'on' => 'S',
		    'or' => 'T',

		    'oid' => 'A',
		    'one' => 'R',
		    'ous' => 'A',

		    'ogen' => 'A',

		    'oidal' => 'A',
		    'oides' => 'A',
		    'otide' => 'A',
		    'ously' => 'A',

		    'oidism' => 'A',

		    'oidally' => 'A',
		    'ousness' => 'A',

		    's' => 'W',

		    "s'" => 'A',

		    'um' => 'U',
		    'us' => 'V',

		    'ward' => 'A',
		    'wise' => 'A',

		    'y' => 'B',

		    'yl' => 'R',

		    'ying' => 'B',
		    'yish' => 'A',

		    "'s" => 'A',
		    );

}
