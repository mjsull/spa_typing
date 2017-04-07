### get_spa_type.py: Get spa types

Version: 0.1.0

License: GPLv3

USAGE: python get_spa_type.py fasta_file.fa

Prints spa type to stdout - egenomics letter combination and then the ridom spa type.

If multiple pcr products are found will print spa types for each product.

Will download sparepeats.fasta and spatypes.txt from the ridom server to repository directory if files not provided or already in directory.
```
optional arguments:

  -r REPEAT_FILE, --repeat_file REPEAT_FILE
                        List of spa repeats
                        (http://spa.ridom.de/dynamic/sparepeats.fasta)
                        
  -o REPEAT_ORDER_FILE, --repeat_order_file REPEAT_ORDER_FILE
                        List spa types and order of repeats
                        (http://spaserver2.ridom.de/dynamic/spatypes.txt)
                        
  -f FASTA, --fasta FASTA
                        Fasta file input.
                        
  --version             show program's version number and exit

```
## How it works

The script searches for 50bp to 5000bp sequences produced by the following primer sets
```
TAAAGACGATCCTTCGGTGAG, CAGCAGTAGTGCCGTTTGCTT
AGACGATCCTTCGGTGAGC, GCTTTTGCAATGTCATTTACTG
ATAGCGTGATTTTGCGGTT, CTAAATATAAATAATGTTGTCACTTGGA
CAACGCAATGGTTTCATCCA, GCTTTTGCAATGTCATTTACTG
```

If an enriched sequence is found by a primer set, subsequent primer sets are not used.

The repeat sequences and repeat orders found on http://spaserver2.ridom.de/ are used to identify the spa type of each enriched sequence.

Ridom spa type and the egenomics repeat sequence are then reported back to the user.



written by mjsull.
