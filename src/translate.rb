#!/usr/bin/env ruby
#
# BioRuby is required.
#
#


def usage
  puts <<HERE
Usage: translate.rb [-options] <RTFINDER output>
  Options:
   -t <n>  : codon table to be used; default = 1 (standard)
   -w <L>  : print windows with length L

HERE
  exit
end


###############################################################################


#
# parse command line options
#
require 'optparse'


$table  = 1
$window = nil
opt = OptionParser.new
opt.on("-w WINDOW") {|w| $window = w.to_i}
opt.on("-t TABLE")  {|t| $table = t.to_i}


begin
  opt.parse!(ARGV)
rescue OptionParser::InvalidOption
  usage
end
usage if ARGV.length == 0


if $window
  if $window % 2 == 0
    STDERR.puts "ERROR: window size should be odd"
    exit
  else
    $window /= 2
  end
end



#
# translate ORFs
#
require 'bio'


ff = Bio::FlatFile.new(Bio::FastaFormat, ARGF)
ff.each_entry do |ent|
  na = Bio::Sequence::NA.new(ent.seq)
  pep = na.translate(1, $table).gsub(/\*/, "U")

  puts ">#{ent.definition}"

  if $window == nil
    puts pep
  else
    ent.definition =~ /stop_codon\=(\d+)/
    pos = $1.to_i - 1

    from = pos - $window
    from = 0 if from < 0

    to = pos + $window
    to = -1 if to >= pep.length

    puts pep[from..to]
  end
end

