# Eimeria_Soay
Metabarcoding of Eimeria spp. in Soay Sheep


## 18S Primers (nuclear genome)

We first test whether 18S primers bind (perfectly, with no mismatches)
to all the Soay Eimeria species' 18S in Aligments. We obtained the
Alignments from Sequences available at NCBI.

Then we test whether the sequence they ampify distinguishes between
the Soay Eimeria species by having at least one nucleotide difference
in the amplified sequence. We list for each of the primers species
that have one difference at least.

## COI Primers (mitochondrial genome)

We align the Sequences from NCBI (only one Soay Eimeria species and
relatives from other ruminants).

I didn't get automatic design to work (because of the high divergence)
and therefore designed primers manually based on the alignemt
"COI_aln.fasta". I store these primers in "COI_Soay_primerMan.txt" in
a format for "https://www.bioinformatics.org/sms2/primer_map.html"
(use the unaligned fasta file "COI_clean.fasta" on the site). See also
"Figures/Primer_Map_COI.png". I also created a fasta file for those
primers: "COI_Soay_Primers.fasta". 


## 28S Primers (nuclear genome)

We align sequences from cow and chicken Eimeria: "28S_aln_clean.fasta"
is the resulting clean alignement file (without gaps in this case). I
"manually" designed primers and they can be found in the primer_map
format in "28S_Soay_primer_map.txt" and fasta in
"28S_Soay_Primers.fasta". A screenshot of the primer map is in
"Figures/Primer_Map_28S.png".


## ORF 470 Primers (apicoplast genome)

We align ORF470 sequences from mouse (!) and chicken Eimeria:
"ORF470_Aln.fasta". I "manually" designed primers and they can be
found in the primer_map format in "ORF470_Soay_primer_map.txt" and
fasta in "ORF470_Soay_Primers.fasta". A screenshot of the primer map
is in "Figures/Primer_Map_ORF470.png".
