# File formats
## summary file

The summary file is ASCII encoded. It consists of two sections—the header and the data. Lines MUST be terminated in the Unix style with a new line (aka line feed) character. Lines in the header section MUST be prefixed with `#`.

The first line of the header section MUST be `# pyhegp summary file version 1`. Subsequent lines of the header section are a list of key-value pairs. Each line MUST be `#`, optional whitespace, the key, a single space character and then the value. The key MUST NOT contain whitespace or control characters, and MUST NOT begin with a `#` character. The value MAY contain whitespace characters, but MUST NOT contain control characters.

The data section is a tab-separated table of numbers. The first line MUST be a header with column labels. Each row corresponds to one SNP. The columns labelled `chromosome`, `position`, `reference`, `mean` and `standard-deviation` contain the chromosome, the position of the SNP on the chromosome, the reference allele, the mean dosage and the standard deviation of the dosage for that SNP. Column labels are case-sensitive.

The `reference` column is optional, and SHOULD be absent in pooled summary files.

Here is an example summary file.
```
# pyhegp summary file version 1
# number-of-samples 100
chromosome	position	reference	mean	standard-deviation
chr11	3200246	G	0.1677	0.3915
chr11	3205355	T	0.1766	0.3599
chr11	3218908	A	0.1818	0.3169
chr11	3414044	A	-0.2420	0.1606
chr11	3460427	C	-0.2408	0.1593
chr11	3462290	A	-0.2383	0.1670
chr11	3462323	A	-0.2388	0.1653
chr11	3462325	T	-0.2386	0.1658
chr11	3462348	T	-0.2371	0.1702
chr11	3464016	A	-0.2400	0.1623
```

## genotype file

The genotype file is a tab-separated values (TSV) file. The first line MUST be a header with column labels. Each row corresponds to one SNP. The columns labelled `chromosome`, `position` and `reference` contain the chromosome, the position on the chromosome and the reference allele for that SNP. Other columns each contain dosage values for one sample. The headers of these columns MUST be their sample identifiers. Column headers are case-sensitive.

the `reference` column is optional, and should be absent in encrypted genotype files.

Here is an example genotype file.
```
chromosome	position	reference	sample1 sample2 sample3 sample4
chr11	3200246	G	0.4043	0.3655	0.3375	-0.614
chr11	3205355	T	0.395	0.3545	0.3325	-0.5421
chr11	3218908	A	0.3977	0.3207	0.3134	-0.4491
chr11	3414044	A	-0.3616	-0.328	-0.3234	0.074
chr11	3460427	C	-0.341	-0.3363	-0.3346	0.073
chr11	3462290	A	-0.3552	-0.3281	-0.3287	0.0912
chr11	3462323	A	-0.3524	-0.3296	-0.3298	0.0874
chr11	3462325	T	-0.3533	-0.3291	-0.3294	0.0885
chr11	3462348	T	-0.361	-0.3244	-0.326	0.0986
chr11	3464016	A	-0.3461	-0.334	-0.3331	0.08
```

## key file

The key file is a tab-separated values (TSV) file with numerical data. There MUST be no column headers.

Here is an example key file.
```
-0.4397501	0.37277363	0.63142739	0.50129352	-0.13290571
0.13568008	-0.13505327	0.2149694	-0.29304762	-0.91173613
0.52289655	0.83464826	0.010733899	-0.17227947	0.012084877
-0.70250368	0.31764876	-0.19654188	-0.60576289	-0.0032338057
-0.14587165	0.21274863	-0.71857058	0.51594477	-0.38848011
```
