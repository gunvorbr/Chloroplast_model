# Chloroplast_model
Model of algal chloroplast developed in MATLAB. The repository also includes scripts developed specifically for working with the chloroplast model.
File iGR774.xml contains a plug-and-play model of the chloroplast of an alga.
The iGR774 model was published in PLoS ONE in 2020 (reference: Røkke, Gunvor Bjerkelund; Hohmann-Marriott, Martin Frank; Almaas, Eivind (2020) An adjustable algal chloroplast plug-and-play model for genome-scale metabolic models. PLoS ONE, vol. 15(2), e0229408, doi: https://doi.org/10.1371/journal.pone.0229408)
The model is tuneable, and can run in three different organism modes - as Nannochloropsis, Chlamydomonas or Phaeodactylum.
Most reactions are common for all three organisms, but certain reactions are organism specific, and these are tagged with either @N_, @C_ or @P_ in front of the reaction ID found in the array model.rxns.
The script changeOrganismMode.m can be used to change organism mode. The default organism mode is Nannochloropsis. changeOrganismMode changes upper and lower (if reversible reaction) bounds of all reactions with an organism tag.
More organism modes can be introduced to the chloroplast model by adding extra reactions to it. New reactions that are added for a specific organism should be tagged in the same manner as described above.
The chloroplast model can be plugged into a larger model structure in need of a chloroplast department. The script plugAndPlay.m will plug the model into an exo-model and change the namespace of the metabolites in the chloroplast model to match the namespace used in the exo-model.
This repository also contains scripts developed specifically for working with the chloroplast model in MATLAB (located in branch Scripts). The scripts all contains descriptions, and can be used alongside the COBRA toolbox.
- Gunvor Røkke, NTNU, 2021
