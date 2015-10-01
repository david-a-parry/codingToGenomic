# CodingToGenomic

CodingToGenomic is a simple program for the conversion of a CDS position to a genomic coordinate using Ensembl identifiers. It is powered by [Ensembl's REST API](http://rest.ensembl.org/).

The following ID types are supported:

* Gene symbols
* Ensembl Gene/Transcript/Protein IDs
* RefSeq IDs
* UniProt IDs
* CCDS IDs
* Entrez Gene IDs

## Installation and Running 

CodingToGenomic is tested on Mac OS X and Linux and should run on most unix-like systems with both Java (v1.7 or higher) and Perl installed. Please consult your OS's documentation for installation instructions for Java (Perl will almost certainly already be installed).

Download the [latest release](https://github.com/gantzgraf/codingToGenomic/releases/latest).

Decompress the downloaded tarball:

    tar -xjvf CodingToGenomic-0.2.tar.bz2

Change directory and run the CodingToGenomic program to get usage information:

    cd CodingToGenomic-0.2
    ./CodingToGenomic

Use the -c/--coordinate option to specify a CDS coordinate to map. Use either the -g/--gene or -t/--transcript options to specify your gene or transcript of interest. The -g/--gene option will attempt to correctly interpret any of the supported ID types (Gene symbols, Ensembl Gene/Transcript/Protein IDs, RefSeq IDs, UniProt IDs, CCDS IDs, Entrez Gene IDs). The -t/--transcript option will always interpret your input as an Ensembl transcript ID, which may be useful for atypical Ensembl transcript IDs such as those for C. elegans.

## Examples

To find the corresponding genomic coordinates for CDS position 100 of ALL identified transcripts for gene TTN (in humans):

    ./CodingToGenomic -g TTN  -c 100

You may also use Ensembl Gene identifiers in place of gene symbols:

    ./CodingToGenomic -g ENSG00000155657  -c 100 

Or other supported identifiers, such as RefSeq IDs:

    ./CodingToGenomic -g NM_001256850.1 -c 100 

To specify a different species other than human use the -s option:

    ./CodingToGenomic -g TTN  -c 100 -s mouse

Use quotes to specify a species name with more than word.

    ./CodingToGenomic -c 100 -g sms-2 -s "Caenorhabditis elegans"

To find the corresponding genomic coordinate for transcript ENST00000540563 CDS position 100:

    ./CodingToGenomic -g ENST00000540563 -c 100

or 

    ./CodingToGenomic -t ENST00000540563 -c 100

If the Ensembl transcript ID is not in the typical "ENS..." format ensure you use the '-t' option to ensure it is interpreted as a transcript ID:

    ./CodingToGenomic -c 100 -t F53H8.4

You may also use other identifiers (CCDS, UniProt, Entrez Gene, RefSeq) with the -g/--gene option to map to the relevant ensembl gene. Other identifiers may also be usable as long as they are not in a format that could be interpretted as one of the supported identifiers and are supported by the Ensembl REST API, so feel free to test different identifiers with the -g/--gene option. 

## Credit

CodingToGenomic was written by David A. Parry (d.a.parry@leeds.ac.uk). It uses [Ensembl's REST API](http://rest.ensembl.org/)

## License

CodingToGenomic is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
