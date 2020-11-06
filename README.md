# Eimeria_Soay
Metabarcoding of Eimeria spp. in Soay Sheep


## 18S Primers

We first test whether 18S primers bind (perfectly, with no mismatches)
to all the Soay Eimeria species' 18S in Aligments. We obtained the
Alignments from Sequences available at NCBI.

Then we test whether the sequence they ampify distinguishes between
the Soay Eimeria species by having at least one nucleotide difference
in the amplified sequence. We list for each of the primers species
that have one difference at least.

## COI Primers

We align the Sequences from NCBI (only one Soay Eimeria species and
relatives from other ruminants).

I didn't get automatic design to work (because of the high divergence)
and therefore designed primers manually based on the alignemt
"COI_aln.fasta". I store these primers in "COI_Soay_primerMan.txt" in
a format for "https://www.bioinformatics.org/sms2/primer_map.html"
(use the unaligned fasta file "COI_clean.fasta" on the site). See also
"Figures/Primer_Map_COI.png". I also created a fasta file for those
primers: "COI_Soay_Primers.fasta". 



