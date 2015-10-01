# CodingToGenomic

CodingToGenomic is a simple program for the conversion of a CDS position to a genomic coordinate using Ensembl identifiers. It is powered by [Ensembl's Rest API](http://rest.ensembl.org/).

## Installation and Running 

CodingToGenomic is tested on Mac OS X and Linux and should run on most unix-like systems with both Java (v1.7 or higher) and Perl installed. Please consult your OS's documentation for installation instructions for Java (Perl will almost certainly already be installed).

Download the [latest release](https://github.com/gantzgraf/codingToGenomic/releases/latest).

Decompress the downloaded tarball:

    tar -xjvf CodingToGenomic-0.1.tar.bz2

Change directory and run the CodingToGenomic program to get usage information:

    cd CodingToGenomic-0.1
    ./CodingToGenomic

## Examples

To find the corresponding genomic coordinate for transcript ENST00000540563 CDS position 100:

    ./CodingToGenomic -t ENST00000540563 -c 100

To find the corresponding genomic coordinates for CDS position 100 of ALL identified transcripts for gene TTN (in humans):

    ./CodingToGenomic -g TTN  -c 100

You may also use Ensembl Gene identifiers in place of gene symbols:

    ./CodingToGenomic -g ENSG00000155657  -c 100 

To specify a different species other than human use the -s option:

    ./CodingToGenomic -g TTN  -c 100 -s mouse

You may also use other identifiers (e.g. CCDS, UniProt, Entrez Gene IDs, RefSeq) with the -g/--gene option to map to the relevant ensembl gene. However, with this feature **all transcripts of the gene will be assessed**, not just the relevant Ensembl transcript, so this feature **SHOULD BE USED WITH CARE**.

## Credit

CodingToGenomic was written by David A. Parry (d.a.parry@leeds.ac.uk). It uses [Ensembl's Rest API](http://rest.ensembl.org/)

## License

CodingToGenomic is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
