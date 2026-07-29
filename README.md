# Norwegian_Ectocarpales_Explororer
The Fehu rune symbol ᚠ means wealth, prosperity, and abundance, a reflection of the dominance of Ectocarpales in terms of brown algal species diversity

The repository hosts a shiny app for exploring Ectocarpales (Phaeophyceae) sequences generated through molecular-morphological surveys conducted in Norway from 2022-2024, funded through the Norwegian Taxonomic Initiative awarded to Kjersti et al.

The app has several features. User inputs include exploring Norwegian survey specific data and/or GenBank data. These data have been curated based on recent literature and taxonomic findings. The User determines which taxa to include by narrowing down the taxonomic fields. A distance model is specified and the tree is computed.

The neighbour joining tree uses coxI data and selected distance model, calculated using the dna.dist function of ape R package. An asterisk indicates a specimen sampled through the Norwegian Taxonomic Initiative (2022-2024). If plot looks messy (adjusting window, data filters), simply hit compute to reset with selected inputs. Note the coxI data are reliable for clustering species, but do not accurately reflect intergeneric relationships. The tree is not bootstrapped or rooted. Some genbank accessions species IDs were updated, see species metadata below for original IDs

A table with specimen information, including previous IDs and geographic locations is also provided. The bioregions are displayed on an interactive map and the user can download the tree as an image or PDF.

This shiny app accompanies the publication of Goerke et al. (2026): An update of the brown algal order Ectocarpales (Phaeophyceae) for Norway, with new records and the reinstatement of Porterinema marina Jaasund. The materials are only uptodate as of the time of publication submission (early 2026).

To reiterate, as noted above, tree performance is not optimized. Users can log concerns through the repo. The shiny app is simply intended to showcase the state of knowledge at the time of publication. We look forward to further advancements in the near future.

We encourage users to adapt this app for other taxonomic studies and to improve user accessibility for molecularly assisted taxonomic surveys. Please see LICENCE.
