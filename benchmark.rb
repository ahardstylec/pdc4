#!/usr/bin/ruby
system("rm bench_erg.txt")
[2, 5, 10].each do |l|
	1..100.times do |t|
		`./QuadCV #{l} #{t} #{t} /home/jakob/Bilder/107_PANA/P1070591.JPG >> bench_erg.txt`
	end
end

f = File.open("bench_erg.txt")
puts "Iterations\tTime\tThreads"
f.each do |l|
	if l.match(/do times (\d+)/) then
		print ($1).to_s.chomp end
	if l.match /time (\d\.\d+)/ then
		print "\t"
		print ($1).to_s.chomp end
	if l.match /There were (\d+)/ then 
		print "\t"
		puts $1.to_s.chomp end
end
