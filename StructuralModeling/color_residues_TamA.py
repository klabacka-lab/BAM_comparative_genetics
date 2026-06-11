# Install the pymol-open-source package to use the pymol package
## conda install pymol-open-source
## the conda-forge chanel should work
from pymol import cmd

cmd.fetch('4C00')
cmd.rotate('y', 80)
cmd.rotate('z', 7)
#cmd.rotate('x', 90)
cmd.select('selection', 'chain A')
cmd.hide('everything', 'selection')
cmd.show('cartoon', 'selection')
cmd.color('black', 'all')
cmd.select('selection', 'chain A and (resi 29 or resi 45 or resi 68 or resi 72 or resi 73 or resi 93 or resi 100 or resi 114 or resi 116 or resi 156 or resi 157 or resi 188 or resi 191 or resi 217 or resi 278 or resi 285 or resi 443 or resi 447 or resi 451 or resi 466 or resi 467 or resi 470 or resi 476 or resi 477 or resi 494 or resi 504 or resi 513 or resi 518 or resi 520 or resi 522 or resi 536 or resi 538 or resi 540 or resi 572 or resi 573)')
cmd.color('red', 'selection')
cmd.show('spheres', 'selection')
cmd.zoom('all', -0.5)
cmd.png('output_image_TamA.png', width=10000, height=10000, dpi=1000)

#cmd.reinitialize()

