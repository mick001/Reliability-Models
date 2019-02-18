# Generates a boolean table of a parallel system

# Generate boolean table
b_table <- gtools::permutations(n=2, r=2, v=c(0,1), repeats.allowed = T)

# Save it
save(b_table, file="bool-table.Rdata")
