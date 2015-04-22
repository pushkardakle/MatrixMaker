# MatrixMaker

*This is a very old code. I have been wanting to refactor it but it has been pending for a long time. If you face any issues/doubts please raise an issue.*

MatrixMaker is a software written in Perl to automate the creation of custom blosum scoring marices for any organism. This approach is useful in the annotation of many hypothetical proteins in organisms with a nucleotide composition deviant from the normal eg. the AT rich E.histolytica/P.falciparum genome

It makes use of UniRef50/90 clusters and the BLOCKS server to create the scoring martices


![Screenshot for MatrixMaker](/home/pushkar/Dropbox/Organizer/Codes/M.Sc/MatrixMaker/Images/matrix_maker_screenshot.png) 

The software accepts the following inputs

* Genus – The genus of the organism
* Species – The species of the organism
* NCBI Taxonomy ID – The NCBI Taxonomy ID for the organism which can be
obtained from - http://www.ncbi.nlm.nih.gov/taxonomy
* Max Per Cluster – This parameter accepts the maximum number of proteins to be
considered per UniRef cluster
* Min Per Cluster – This parameter accepts the minimum number of proteins that a
UniRef cluster should have to be considered
* Total Number of Clusters – This is the total number of UniRef clusters that will
be considered for the creation of the matrix. A lower number would increase
speed but compromise accuracy.
* Taxons to Select – This parameter accepts a list of scientific names separated by
newline for the organisms to which the matrix will be limited. If kept blank there
will be no taxonomic limitations for the creation of the matrix.
* Proxy – If enabled the network connections will be made through specified proxy.
* Prioritization – This parameter if enables emulates the selection of proteins from
taxonomically diverse organisms for the creation of the matrix.
* UniRefNo – Gives the option of UniRef50/UniRef90 clusters to be chosen for the
creation of the matrix.



