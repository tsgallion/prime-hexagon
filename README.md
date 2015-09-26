# prime-hexagon
Prime numbers and their positions in the prime hexagon

This code plots the positons in the central triangle of the prime hexagon. 

What is the prime hexagon, you may ask? TODO...

##  Computing a prime number's position on the prime hexagon

Initial prime within the triangle is at position 2 and is prime 5.  Numbers at 2,4,0
are dependent on status of those at 1,3, and 5 for spin determination for numbers
immediately following.  (Because they are n type primes and defined as residing
in the last hexagon.  The following number, n+1, spins with the n number, unless
it is followed by another prime, then spins with the second prime.
A sketch of the inner triangle:


           <- 3 ->
           
           2      4
           
       <-1    0     5->
         |          |
         v          v
         
