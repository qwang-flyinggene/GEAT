# GEAT 

##What is it?
GEAT is a java-based **G**enomic **E**vents **A**nalysis **T**ool aimed to performe analysis, visualization and modelling of genomic events data obtained from high throughput platforms.

##Download & Install

```$ git clone "https://github.com/geatools/geat.git/"```

##Requirements / Dependencies##
- **Java**    
- **BLAST+**     
> *GEAT employs blastn to performe two short seq alignment. To setsup BLAST+, Go to the [BLAS T download page](http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)*.
- **Perl**    
> *GEAT currently performes seq QC check based on included [Prinseq](http://prinseq.sourceforge.net/) which is developed by [Perl](https://www.perl.org/).*

## Execution/Command

java -Xmx16G -jar GEAT.jar [options]

