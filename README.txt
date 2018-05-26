Once upon a time my dad worked on supercomputers at Los Alamos National Laboratory and sometimes he would spend a minute running a little test to see about how fast a new supercomputer was. That test is here as flops.c very much as I received it in 1992. Your phone is probably much faster than the fastest supercomputer of 1992.
This isn't a very good benchmark, but it's simple and easy. It runs a few basic numeric algorithms with a known number of adds, subtracts, multiplies, and divides, and figures out how many floating point operations per second your CPU can do. It doesn't take advantage of multi-core or SIMD instructions. It doesn't exercise the memory system and probably all fits in L1 cache. It just tests how fast your CPU can do math (and that your compiler isn't terrible at making that happen).
Over the years I transliterated flops.c into other languages to test them and their compilers and interpreters. Python has a pretty slow interpreter, around 1-5% of the speed of C. Javascript got amazingly good and can run 80-90% of the speed of C. Java got to that speed around 2007 or 2007. The Go compiler is surprisingly good for a newer language. Julia is a newer language that should have the potential for running full speed but apparently still needs some tweaking (as of 2018-05).


Running


cc flops.c && ./a.out
javac flops.java && java flops
python flops.py
julia flops.jl
go run flops.go
open flops.html in your browser to run flops.js
