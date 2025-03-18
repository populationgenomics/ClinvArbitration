These files were downloaded from the Ensembl FTP sites 6th February 2025. They are stored off-chain using git lfs, and I
planned to make that the default distribution method for the files. However, I've since realised that git-lfs has quite
a limited default quota, and you sourcing these files separately from ensembl is far simpler.

To download these files, either install git-lfs locally (which will automatically populate them when this repository is
cloned, assuming we are within data caps), or visit the Ensembl FTP sites:

GRCh38:

https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/

or for GRCh37

https://ftp.ensembl.org/pub/grch37/release-113/gff3/homo_sapiens/

Once downloaded, rename them to a standardised <ASSEMBLY>.gff3.gz name format:

- GRCh38.gff3.gz
- GRCh37.gff3.gz

The committed versions of these files were originally sourced from Ensembl data release 113, and have had their names truncated:

- GRCh38.gff3.gz was previously Homo_sapiens.GRCh38.113.gff3.gz
- GRCh37.gff3.gz was previously Homo_sapiens.GRCh37.87.gff3.gz
