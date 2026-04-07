reverse_complement_vector <-
function(dna_seqs) {
  DNAStringSet(dna_seqs) |> reverseComplement() |> as.character()
}
