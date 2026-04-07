reverse_complement <-
function(dna_seq) {
  sapply(dna_seq, function(x) {
    DNAString(x) |> reverseComplement() |> as.character()
  })
}
