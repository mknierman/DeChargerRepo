![deCharger](https://github.com/mknierman/mknierman.github.io/blob/gh-pages/Decharger_icon-small.png)
# deCharger
Mass Spectrometry Decharged Viewer
Developed originally by Tony Major, Michael Knierman, and Anvesh Kodumuri

Requirements: Windows PC (64-bit) with .NET Framework 4.8 or higher

Installation
Download the zip file and run the decharger_setup.exe in the root directory.

deCharger was originally intended as a proof of Concept in 2012.  The intent was to make full use of the mass accuracy and isotopic resolution available in MS2 spectra from newer instrumentation.  For the higher mass ions found in many top down experiments, the charge detection of contemporary tools frequently failed to correctly identify the charge state and frequently failed to identify the monoisotopic mass.  DeCharger was an attempt to address this.  As time went on, we added more features to the point where it can find de novo sequence tags, and search those tags in a database and return the protein identification results.  The tool was originally coded for CID only data, but because of the nature of how it works, it works to some extent with other activations.  However, a refactor of this is needed to properly handle c/z ion fragments. 

To corectly handle c /z ions in the validation workbench add [17] to the N terminal and [-16] to the C terminal in the sequence.

deCharger can handle basic oligonucleotides (RNA, DNA bases) and can somewhat handle oligo ms/ms spectra.  Set the isotope model in the config page (gear icon) to use oligo "averagine" for the isotope model.  In the match list tab select DNA or RNA as the base mass list.  You can add to the mass list by adding a symbol in the format of a capitol letter followed by a small case letter [Aa] <tab> the  4 decimal monoisotopic mass of the base.  The headers of the fragment ion table in the validation tab are still for peptides.  To make the mass match in the validation tab,  subtract 80 from the 5' end (N-term peptide).  The y ions are w ions.   To make the b ions into A-B ions add 98.02 to the 5' end
