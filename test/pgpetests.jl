
@test modeorthonorm("Hermite",100,1e-10) #test orthonormality
@test modeorthonorm("Hermite",100,1e-10,3) #test scaling with frequency
@test modeorthonorm("Hermite",370,1e-8) #test maximum cutoff
