# Install the pymol-open-source package to use the pymol package
## conda install pymol-open-source
## the conda-forge chanel should work
from pymol import cmd

cmd.fetch('8ADG')
cmd.rotate('y', 65)
cmd.rotate('z', 45)
#cmd.rotate('x', 45)
cmd.select('selection', 'chain B or chain C or chain D or chain E or chain F')
cmd.hide('everything', 'selection')
cmd.color('black', 'all')
cmd.select('selection', 'chain A and (resi 50 or resi 183 or resi 233 or resi 315 or resi 355 or resi 367 or resi 419 or resi 447 or resi 450 or resi 589 or resi 646 or resi 655 or resi 660 or resi 661 or resi 737 or resi 739 or resi 770 or resi 806)')
cmd.color('red', 'selection')
cmd.show('spheres', 'selection')
cmd.zoom('all', -0.5)
cmd.png('output_image.png', width=4000, height=4000, dpi=2000)

#cmd.reinitialize()

