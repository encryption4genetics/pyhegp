# File formats
## summary file

The summary file is ASCII encoded. It consists of two sections—the header and the data. Lines MUST be terminated in the Unix style with a new line (aka line feed) character. Lines in the header section MUST be prefixed with `#`.

The first line of the header section MUST be `# pyhegp summary file version 1`. Subsequent lines of the header section are a list of key-value pairs. Each line MUST be `#`, optional whitespace, the key, a single space character and then the value. The key MUST NOT contain whitespace or control characters, and MUST NOT begin with a `#` character. The value MAY contain whitespace characters, but MUST NOT contain control characters.

The data section is a space separated table of numbers. The first line of the data section is a vector of means—one for each SNP. The second line is a vector of standard deviations—one for each SNP.

Here is an example summary file.
`TODO: Add example.`
