# Auto-generated PyMOL script for BioViz Mutation Analysis
cmd.hide('everything', 'target_protein')
cmd.show('cartoon', 'target_protein')
cmd.color('grey80', 'target_protein')

# Select and highlight mutated residues
cmd.select('mutated_sites', 'resi 1+2+3+10+17 and target_protein')
cmd.show('sticks', 'mutated_sites')
cmd.color('firebrick', 'mutated_sites')

# Set viewport dimensions, focus camera framing, and render image
cmd.viewport(1000, 1000)
cmd.zoom('mutated_sites', buffer=5.0)
cmd.ray(1000, 1000)
cmd.png('C:/Users/jyam4/PycharmProjects/GenomePilot/src/bioviz/mutations/PyMOL_rendering/mutation_render.png')