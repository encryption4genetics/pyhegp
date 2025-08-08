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
`TODO: Add example.`

## key file

The key file is a tab-separated values (TSV) file with numerical data. There MUST be no column headers.

Here is an example key file.
`TODO: Add example.`
